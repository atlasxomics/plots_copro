new_data_signal()

# The heatmap cell defines the shared selector and neighborhood helpers before
# exiting when this view is active.
_neigh_display_value = (
    neigh_display.value if "neigh_display" in globals() else "heatmap"
)
if _neigh_display_value != "radial":
    exit()

w_text_output(
    content="""
## Radial Neighborhood Enrichment

Each subplot is one focal cluster. Spokes point to neighboring clusters; every
spoke reaches the perimeter, **red/blue** encodes positive/negative enrichment,
and thickness encodes the absolute z-score. Hover a spoke for its value.
"""
)

if adata_ge is None:
    w_text_output(
        content="No ATAC gene-accessibility data selected...",
        appearance={"message_box": "warning"},
    )
    exit()


def build_radial_neighborhood_fig(
    z_matrix,
    labels,
    title,
    threshold=0.0,
    show_self=False,
    shared_zmax=None,
):
    """Build a grid of fixed-radius polar spoke plots."""
    n_clusters = len(labels)
    if shared_zmax is None:
        off_diagonal = np.abs(z_matrix).copy()
        np.fill_diagonal(off_diagonal, 0.0)
        finite = off_diagonal[np.isfinite(off_diagonal)]
        shared_zmax = float(np.nanmax(finite)) if finite.size else 1.0
    if not np.isfinite(shared_zmax) or shared_zmax <= 0:
        shared_zmax = 1.0

    width_min, width_max = 0.04, 0.55
    ncols = min(4, max(1, n_clusters))
    nrows = math.ceil(n_clusters / ncols)
    fig = make_subplots(
        rows=nrows,
        cols=ncols,
        specs=[[{"type": "polar"} for _ in range(ncols)] for _ in range(nrows)],
        subplot_titles=labels,
        horizontal_spacing=0.04,
        vertical_spacing=0.11,
    )

    for focal_index, focal in enumerate(labels):
        radii = []
        angles = []
        colors = []
        widths = []
        zscores = []
        for neighbor_index, neighbor in enumerate(labels):
            if not show_self and focal_index == neighbor_index:
                continue
            zscore = float(z_matrix[focal_index, neighbor_index])
            if not np.isfinite(zscore) or abs(zscore) < threshold:
                continue
            denominator = shared_zmax - threshold
            fraction = (
                (abs(zscore) - threshold) / denominator if denominator > 0 else 1.0
            )
            fraction = min(1.0, max(0.0, fraction))
            radii.append(1.0)
            angles.append(str(neighbor))
            colors.append("#C33530" if zscore >= 0 else "#282E66")
            widths.append(width_min + fraction * (width_max - width_min))
            zscores.append(zscore)

        fig.add_trace(
            go.Barpolar(
                theta=angles,
                r=radii,
                marker=dict(color=colors, line=dict(width=0)),
                width=widths,
                customdata=zscores,
                showlegend=False,
                hovertemplate=(
                    f"focal {focal}<br>neighbor %{{theta}}"
                    "<br>z-score %{customdata:.2f}<extra></extra>"
                ),
            ),
            row=focal_index // ncols + 1,
            col=focal_index % ncols + 1,
        )

    fig.update_polars(
        radialaxis=dict(range=[0, 1], showticklabels=False, ticks="", showline=False),
        angularaxis=dict(
            direction="clockwise",
            rotation=90,
            categoryorder="array",
            categoryarray=[str(label) for label in labels],
            tickfont=dict(size=7),
        ),
        bgcolor="white",
    )
    fig.update_layout(
        title=dict(text=title, x=0.5, xanchor="center", font=dict(size=16)),
        width=250 * ncols,
        height=270 * nrows,
        margin=dict(l=20, r=20, t=80, b=20),
    )
    for annotation in fig.layout.annotations:
        annotation.font.size = 12
        annotation.yshift = 10
    return fig


selected_cluster_key = neigh_cluster_key.value
radial_precomputed = load_precomputed_neighborhood_groups(
    adata_ge, selected_cluster_key
)
excluded_group_keys = set(available_neighborhood_clusters) | set(NA_KEYS)
radial_group_keys = [
    key for key in get_groupable_obs_keys(adata_ge) if key not in excluded_group_keys
]
for key in radial_precomputed:
    if key in adata_ge.obs and key not in radial_group_keys:
        radial_group_keys.insert(0, key)

radial_group_by = w_select(
    label="subplot groups",
    key="radial_group_by",
    default="all",
    options=tuple(["all"] + radial_group_keys),
    appearance={"help_text": "Facet radial plots by a categorical observation."},
)
radial_mode = w_select(
    label="value metric",
    key="radial_mode",
    default="zscore",
    options=("zscore", "count"),
    appearance={
        "help_text": (
            "Metric included in the table; spoke color and thickness always use z-score."
        )
    },
)
radial_threshold = w_text_input(
    label="z-score threshold",
    key="radial_threshold",
    default="0",
    appearance={
        "help_text": "Only spokes at or above this absolute z-score are drawn."
    },
)
radial_show_self = w_checkbox(
    label="Show self-enrichment spoke",
    key="radial_show_self",
    default=False,
)
w_row(
    items=[
        radial_group_by,
        radial_mode,
        radial_threshold,
        radial_show_self,
    ]
)

if radial_group_by.value not in (None, "all", *radial_precomputed.keys()):
    w_text_output(
        content=(
            "This annotation was added after the workflow ran. Spatial "
            "neighborhoods will be computed when first displayed."
        ),
        appearance={"message_box": "warning"},
    )

radial_button = w_button(label="Update Radial Plots", key="radial_button")

if radial_group_by.value is not None and radial_button.value:
    try:
        threshold = (
            float(radial_threshold.value)
            if radial_threshold.value not in (None, "")
            else 0.0
        )
    except (TypeError, ValueError):
        threshold = 0.0
        w_text_output(
            content="z-score threshold must be numeric; defaulting to 0.",
            appearance={"message_box": "warning"},
        )
    threshold = max(0.0, threshold)
    show_self = bool(radial_show_self.value)
    mode = radial_mode.value
    result_key = neighborhood_result_key(selected_cluster_key)

    if radial_group_by.value == "all":
        if result_key in adata_ge.uns:
            plot_adatas = {"all": adata_ge}
            status = "Using workflow-precomputed neighborhood enrichment."
        else:
            if selected_cluster_key not in neighborhood_all_results:
                w_text_output(
                    content="Computing neighborhoods for all spots...",
                    appearance={"message_box": "info"},
                )
                submit_widget_state()
                lightweight = make_lightweight_neighborhood_adata(
                    adata_ge, selected_cluster_key
                )
                compute_neighborhood_enrichment(lightweight, selected_cluster_key)
                neighborhood_all_results[selected_cluster_key] = lightweight
            plot_adatas = {"all": neighborhood_all_results[selected_cluster_key]}
            status = "Using locally computed neighborhood enrichment."
        w_text_output(content=status, appearance={"message_box": "info"})
        submit_widget_state()
        group_label = "all"
    else:
        group_label = radial_group_by.value
        cache_key = (selected_cluster_key, group_label)
        if group_label in radial_precomputed:
            neighborhood_filtered_groups[cache_key] = radial_precomputed[group_label]
        elif cache_key not in neighborhood_filtered_groups:
            neighborhood_filtered_groups[cache_key] = {}
        plot_adatas = neighborhood_filtered_groups[cache_key]

        subgroup_values = [
            str(value) for value in adata_ge.obs[group_label].dropna().unique()
        ]
        for subgroup in subgroup_values:
            if subgroup not in plot_adatas:
                plot_adatas[subgroup] = make_lightweight_neighborhood_adata(
                    adata_ge, selected_cluster_key, group_label, subgroup
                )
            if result_key not in plot_adatas[subgroup].uns:
                w_text_output(
                    content=f"Computing spatial neighborhoods for {subgroup}...",
                    appearance={"message_box": "info"},
                )
                submit_widget_state()
                compute_neighborhood_enrichment(
                    plot_adatas[subgroup], selected_cluster_key
                )

    matrices = {}
    shared_zmax = 0.0
    for subgroup, subgroup_adata in plot_adatas.items():
        z_matrix, labels = ordered_neighborhood_matrix(
            subgroup_adata, selected_cluster_key, "zscore"
        )
        mode_matrix, mode_labels = ordered_neighborhood_matrix(
            subgroup_adata, selected_cluster_key, mode
        )
        if labels != mode_labels:
            raise ValueError("Neighborhood metric labels do not align.")
        matrices[str(subgroup)] = (z_matrix, mode_matrix, labels)
        off_diagonal = np.abs(z_matrix).copy()
        np.fill_diagonal(off_diagonal, 0.0)
        finite = off_diagonal[np.isfinite(off_diagonal)]
        if finite.size:
            shared_zmax = max(shared_zmax, float(np.nanmax(finite)))
    if shared_zmax <= 0:
        shared_zmax = 1.0

    radial_figures = []
    radial_rows = []
    for subgroup in sort_group_categories(list(matrices)):
        z_matrix, mode_matrix, labels = matrices[subgroup]
        title = (
            f"All spots: {selected_cluster_key} Radial Neighborhood Enrichment"
            if group_label == "all"
            else (
                f"{group_label} = {subgroup}: {selected_cluster_key} "
                "Radial Neighborhood Enrichment"
            )
        )
        radial_figures.append(
            build_radial_neighborhood_fig(
                z_matrix,
                labels,
                title,
                threshold=threshold,
                show_self=show_self,
                shared_zmax=shared_zmax,
            )
        )
        for focal_index, focal in enumerate(labels):
            for neighbor_index, neighbor in enumerate(labels):
                radial_rows.append(
                    {
                        "subgroup": subgroup,
                        "focal_cluster": focal,
                        "neighbor_cluster": neighbor,
                        "zscore": z_matrix[focal_index, neighbor_index],
                        mode: mode_matrix[focal_index, neighbor_index],
                    }
                )

    for index, figure in enumerate(radial_figures):
        w_plot(source=figure, key=f"radial_plot_{index}")
    radial_enrichment_data = pd.DataFrame(radial_rows)
    w_table(
        label="Radial neighborhood enrichment values",
        source=radial_enrichment_data,
        key="radial_enrichment_data",
    )
