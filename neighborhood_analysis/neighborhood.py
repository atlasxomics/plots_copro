new_data_signal()

import plotly.graph_objects as go
from plotly.subplots import make_subplots


w_text_output(
    content="""
# Neighborhood Analysis

Explore spatial neighborhood enrichment for the **CoPro**, **RNA**, or **ATAC**
cluster assignments. Results precomputed by the SpatialGlue workflow are used
when available; custom annotations are computed when first displayed.

<details>
<summary><i>details</i></summary>

Each heatmap shows how often spots from one cluster neighbor spots from another
compared with chance. **zscore** is the standardized enrichment and is best for
comparisons; **count** is the raw neighborhood count. Plots can include all
spots or be split by a categorical observation such as sample or condition.

Use **display type** to switch between this matrix view and the radial spoke
view below. Optional color-scale limits apply to every heatmap in the view.
</details>
"""
)

if adata_ge is None:
    w_text_output(
        content="No ATAC gene-accessibility data selected...",
        appearance={"message_box": "warning"},
    )
    exit()


NHOOD_SUFFIX = "_nhood_enrichment"
GROUPED_NHOOD_SUFFIX = "_nhood_enrichment_by_group"


def neighborhood_result_key(cluster_key):
    return f"{cluster_key}{NHOOD_SUFFIX}"


def grouped_neighborhood_result_key(cluster_key):
    return f"{cluster_key}{GROUPED_NHOOD_SUFFIX}"


def neighborhood_cluster_keys(adata):
    """Return cluster annotations, preferring workflow-precomputed results."""
    precomputed = []
    for key in adata.uns:
        if key.endswith(NHOOD_SUFFIX) and not key.endswith(GROUPED_NHOOD_SUFFIX):
            cluster_key = key[: -len(NHOOD_SUFFIX)]
            if cluster_key in adata.obs:
                precomputed.append(cluster_key)

    preferred = ["CoPro_cluster", "RNA_cluster", "ATAC_cluster"]
    ordered = [key for key in preferred if key in precomputed]
    ordered.extend(sorted(key for key in precomputed if key not in ordered))

    for key in get_cluster_keys(adata):
        if key not in ordered:
            ordered.append(key)
    return ordered


def load_precomputed_neighborhood_groups(adata, cluster_key):
    """Rebuild lightweight AnnData objects from the workflow's compact schema."""
    root = adata.uns.get(grouped_neighborhood_result_key(cluster_key))
    if not isinstance(root, dict) or int(root.get("schema_version", -1)) != 1:
        return {}

    result = {}
    for group_entry in root.get("groups", {}).values():
        group_key = str(group_entry["group_key"])
        subgroups = {}
        for subgroup_entry in group_entry.get("subgroups", {}).values():
            group_value = str(subgroup_entry["group_value"])
            categories = pd.Index(
                np.asarray(subgroup_entry["cluster_categories"]).astype(str)
            )
            obs = pd.DataFrame(
                {
                    cluster_key: pd.Categorical(
                        categories,
                        categories=categories,
                        ordered=True,
                    )
                }
            )
            neighborhood_adata = anndata.AnnData(obs=obs)
            neighborhood_adata.uns[neighborhood_result_key(cluster_key)] = {
                "zscore": np.asarray(subgroup_entry["zscore"]),
                "count": np.asarray(subgroup_entry["count"]),
            }
            subgroups[group_value] = neighborhood_adata
        result[group_key] = subgroups
    return result


def make_lightweight_neighborhood_adata(adata, cluster_key, group=None, subgroup=None):
    """Copy only annotations and coordinates, never the feature matrix."""
    if DEFAULT_DATASET_OBSM_KEY not in adata.obsm:
        raise KeyError(
            f"Expected `{DEFAULT_DATASET_OBSM_KEY}` in adata.obsm "
            "(created while loading data)."
        )

    mask = np.ones(adata.n_obs, dtype=bool)
    if group is not None:
        mask &= (adata.obs[group].astype(str) == str(subgroup)).to_numpy()
    mask &= adata.obs[cluster_key].notna().to_numpy()

    obs_keys = [cluster_key]
    if "sample" in adata.obs.columns and "sample" not in obs_keys:
        obs_keys.append("sample")
    subset_obs = adata.obs.loc[mask, obs_keys].copy()
    for obs_key in obs_keys:
        if isinstance(subset_obs[obs_key].dtype, pd.CategoricalDtype):
            subset_obs[obs_key] = subset_obs[obs_key].cat.remove_unused_categories()

    lightweight = anndata.AnnData(obs=subset_obs)
    lightweight.obsm[DEFAULT_DATASET_OBSM_KEY] = np.asarray(
        adata.obsm[DEFAULT_DATASET_OBSM_KEY]
    )[mask].copy()
    return lightweight


def compute_neighborhood_enrichment(adata, cluster_key):
    """Compute Squidpy enrichment on a lightweight AnnData object."""
    try:
        from squidpy.gr import nhood_enrichment, spatial_neighbors
    except ImportError as error:
        raise RuntimeError(
            "Squidpy is required to compute neighborhoods for annotations that "
            "were not precomputed by the workflow."
        ) from error

    if not isinstance(adata.obs[cluster_key].dtype, pd.CategoricalDtype):
        adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")
    else:
        adata.obs[cluster_key] = adata.obs[cluster_key].cat.remove_unused_categories()

    sample_key = "sample" if "sample" in adata.obs else None
    if sample_key is not None:
        if not isinstance(adata.obs[sample_key].dtype, pd.CategoricalDtype):
            adata.obs[sample_key] = adata.obs[sample_key].astype("category")
        else:
            adata.obs[sample_key] = adata.obs[sample_key].cat.remove_unused_categories()

    result_key = neighborhood_result_key(cluster_key)
    if adata.n_obs < 2 or len(adata.obs[cluster_key].cat.categories) < 2:
        n_clusters = len(adata.obs[cluster_key].cat.categories)
        zeros = np.zeros((n_clusters, n_clusters), dtype=float)
        adata.uns[result_key] = {"zscore": zeros, "count": zeros.copy()}
        return adata

    spatial_neighbors(
        adata,
        spatial_key=DEFAULT_DATASET_OBSM_KEY,
        coord_type="grid",
        n_neighs=4,
        n_rings=1,
        library_key=sample_key,
    )
    nhood_enrichment(
        adata,
        cluster_key=cluster_key,
        library_key=sample_key,
        seed=42,
    )
    return adata


def ordered_neighborhood_matrix(adata, cluster_key, mode):
    """Return a matrix and consistently ordered string cluster labels."""
    result_key = neighborhood_result_key(cluster_key)
    if result_key not in adata.uns:
        raise KeyError(f"Neighborhood result `{result_key}` was not found.")

    series = adata.obs[cluster_key]
    if not isinstance(series.dtype, pd.CategoricalDtype):
        series = series.astype("category")
    categories = [str(value) for value in series.cat.categories]
    matrix = np.asarray(adata.uns[result_key][mode], dtype=float)
    if matrix.shape != (len(categories), len(categories)):
        raise ValueError(
            f"Neighborhood matrix {result_key!r} has shape {matrix.shape}, but "
            f"{cluster_key!r} has {len(categories)} categories."
        )

    labels = sort_group_categories(categories)
    order = [categories.index(label) for label in labels]
    return matrix[np.ix_(order, order)], labels


def neighborhood_colorscale(vmin, vmax):
    if vmin >= 0:
        return [[0, "white"], [1, "red"]]
    if vmax <= 0:
        return [[0, "blue"], [1, "white"]]
    midpoint = abs(vmin) / (abs(vmin) + abs(vmax))
    return [[0, "blue"], [midpoint, "white"], [1, "red"]]


def plot_neighborhood_heatmaps(
    group_adatas, cluster_key, title, mode, vmin=None, vmax=None
):
    """Build one or more heatmaps with a shared scale and a labeled data table."""
    group_names = sort_group_categories([str(name) for name in group_adatas])
    matrices = {}
    labels_by_group = {}
    for group_name in group_names:
        matrix, labels = ordered_neighborhood_matrix(
            group_adatas[group_name], cluster_key, mode
        )
        matrices[group_name] = matrix
        labels_by_group[group_name] = labels

    finite_mins = [np.nanmin(value) for value in matrices.values() if value.size]
    finite_maxs = [np.nanmax(value) for value in matrices.values() if value.size]
    if vmin is None:
        vmin = float(min(finite_mins)) if finite_mins else 0.0
    if vmax is None:
        vmax = float(max(finite_maxs)) if finite_maxs else 0.0
    colorscale = neighborhood_colorscale(vmin, vmax)

    n_groups = len(group_names)
    ncols = min(4, max(1, n_groups))
    nrows = math.ceil(n_groups / ncols)
    cell_width = 520 if n_groups <= 4 else 380
    cell_height = 480 if n_groups <= 4 else 350
    fig = make_subplots(
        rows=nrows,
        cols=ncols,
        subplot_titles=group_names,
        horizontal_spacing=0.06,
        vertical_spacing=0.08,
    )

    table_rows = []
    for index, group_name in enumerate(group_names):
        row = index // ncols + 1
        col = index % ncols + 1
        matrix = matrices[group_name]
        labels = labels_by_group[group_name]
        fig.add_trace(
            go.Heatmap(
                z=matrix,
                x=labels,
                y=labels,
                colorscale=colorscale,
                zmin=vmin,
                zmax=vmax,
                showscale=index == 0,
                hoverongaps=False,
                colorbar=dict(title=mode, tickformat=".2f") if index == 0 else None,
                hovertemplate=(
                    "focal %{y}<br>neighbor %{x}<br>"
                    + mode
                    + " %{z:.2f}<extra></extra>"
                ),
            ),
            row=row,
            col=col,
        )
        fig.update_yaxes(autorange="reversed", row=row, col=col)
        for focal_index, focal in enumerate(labels):
            for neighbor_index, neighbor in enumerate(labels):
                table_rows.append(
                    {
                        "subgroup": group_name,
                        "focal_cluster": focal,
                        "neighbor_cluster": neighbor,
                        mode: matrix[focal_index, neighbor_index],
                    }
                )

    fig.update_layout(
        title=dict(text=title, x=0.5, xanchor="center"),
        width=cell_width * ncols,
        height=cell_height * nrows,
        margin=dict(l=70, r=90, t=90, b=50),
        plot_bgcolor="white",
    )
    return fig, pd.DataFrame(table_rows)


def parse_optional_float(value, label):
    if value is None or str(value).strip() == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        w_text_output(
            content=f"{label} must be numeric; using automatic scaling.",
            appearance={"message_box": "warning"},
        )
        return None


available_neighborhood_clusters = neighborhood_cluster_keys(adata_ge)
if not available_neighborhood_clusters:
    w_text_output(
        content="No cluster annotations are available for neighborhood analysis.",
        appearance={"message_box": "warning"},
    )
    exit()

default_neighborhood_cluster = choose_default_option(
    available_neighborhood_clusters,
    preferred="CoPro_cluster",
)
neigh_display = w_select(
    label="display type",
    key="neigh_display",
    default="heatmap",
    options=("heatmap", "radial"),
    appearance={"help_text": "Switch between matrix and radial spoke views."},
)
neigh_cluster_key = w_select(
    label="cluster annotation",
    key="neigh_cluster_key",
    default=default_neighborhood_cluster,
    options=tuple(available_neighborhood_clusters),
    appearance={"help_text": "Choose the clustering whose neighborhoods are compared."},
)
w_row(items=[neigh_display, neigh_cluster_key])

if neigh_display.value != "heatmap":
    exit()

selected_cluster_key = neigh_cluster_key.value
precomputed_neighborhoods = load_precomputed_neighborhood_groups(
    adata_ge, selected_cluster_key
)
excluded_group_keys = set(available_neighborhood_clusters) | set(NA_KEYS)
neighbor_groups = [
    key for key in get_groupable_obs_keys(adata_ge) if key not in excluded_group_keys
]
for key in precomputed_neighborhoods:
    if key in adata_ge.obs and key not in neighbor_groups:
        neighbor_groups.insert(0, key)

neigh_group_by = w_select(
    label="subplot groups",
    key="neigh_group_by",
    default="all",
    options=tuple(["all"] + neighbor_groups),
    appearance={"help_text": "Facet by a categorical observation."},
)
neigh_mode = w_select(
    label="displayed data",
    key="neigh_mode",
    default="zscore",
    options=("zscore", "count"),
)
neigh_scale_max = w_text_input(
    label="colorscale maximum", key="neigh_scale_max", default=None
)
neigh_scale_min = w_text_input(
    label="colorscale minimum", key="neigh_scale_min", default=None
)
w_row(items=[neigh_group_by, neigh_mode, neigh_scale_max, neigh_scale_min])

if neigh_group_by.value not in (None, "all", *precomputed_neighborhoods.keys()):
    w_text_output(
        content=(
            "This annotation was added after the workflow ran. Spatial "
            "neighborhoods will be computed when first displayed."
        ),
        appearance={"message_box": "warning"},
    )

neigh_button = w_button(label="Update Neighborhood Plots", key="neigh_button")

if neigh_group_by.value is not None and neigh_button.value:
    vmin = parse_optional_float(neigh_scale_min.value, "Colorscale minimum")
    vmax = parse_optional_float(neigh_scale_max.value, "Colorscale maximum")
    result_key = neighborhood_result_key(selected_cluster_key)

    if neigh_group_by.value == "all":
        if result_key in adata_ge.uns:
            plot_adatas = {"all": adata_ge}
            status = "Using workflow-precomputed neighborhood enrichment."
        else:
            if selected_cluster_key not in neighborhood_all_results:
                status = "Computing neighborhoods for all spots..."
                w_text_output(content=status, appearance={"message_box": "info"})
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
        title = f"All spots: {selected_cluster_key} Neighborhood Enrichment"
    else:
        group = neigh_group_by.value
        cache_key = (selected_cluster_key, group)
        if group in precomputed_neighborhoods:
            neighborhood_filtered_groups[cache_key] = precomputed_neighborhoods[group]
        elif cache_key not in neighborhood_filtered_groups:
            neighborhood_filtered_groups[cache_key] = {}

        plot_adatas = neighborhood_filtered_groups[cache_key]
        subgroup_values = [
            str(value) for value in adata_ge.obs[group].dropna().unique()
        ]
        for subgroup in subgroup_values:
            if subgroup not in plot_adatas:
                plot_adatas[subgroup] = make_lightweight_neighborhood_adata(
                    adata_ge, selected_cluster_key, group, subgroup
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
        title = f"{selected_cluster_key} Neighborhoods by {group}"

    neigh_heatmap, neigh_data = plot_neighborhood_heatmaps(
        plot_adatas,
        selected_cluster_key,
        title,
        neigh_mode.value,
        vmin=vmin,
        vmax=vmax,
    )
    w_plot(source=neigh_heatmap, key="neigh_heatmap")
    w_table(source=neigh_data, key="neigh_data")
