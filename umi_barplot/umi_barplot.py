w_text_output(content="""
## Mean Counts Per Spot

Generate side-by-side WT and GE barplots of mean counts per spot grouped by
SpatialGlue cluster for a selected feature.
""")

new_data_signal()
if adata_rna is None and adata_ge is None:
    w_text_output(
        content="No data loaded...",
        appearance={"message_box": "warning"},
    )
    fig = None
    submit_widget_state()
    exit()

count_objects = []
if adata_rna is not None:
    count_objects.append({
        "label": "WT",
        "adata": adata_rna,
        "features": list(available_genes),
        "object_name": rna_object_name,
    })
if adata_ge is not None:
    count_objects.append({
        "label": "GE",
        "adata": adata_ge,
        "features": list(available_ge_features),
        "object_name": ge_object_name,
    })

for obj in count_objects:
    if not obj["features"]:
        w_text_output(
            content=f"No features were found in `{obj['object_name']}.h5ad`.",
            appearance={"message_box": "warning"},
        )
        submit_widget_state()
        exit()

if len(count_objects) > 1:
    feature_sets = [set(obj["features"]) for obj in count_objects]
    shared_features = set.intersection(*feature_sets)
    feature_options = [
        feature
        for feature in count_objects[0]["features"]
        if feature in shared_features
    ]
else:
    feature_options = count_objects[0]["features"]

if not feature_options:
    w_text_output(
        content="No shared WT/GE features were found for side-by-side plotting.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

cluster_options_by_object = []
for obj in count_objects:
    cluster_options = get_cluster_keys(obj["adata"])
    if not cluster_options:
        w_text_output(
            content=f"No cluster-like metadata columns were found in `{obj['object_name']}.h5ad`.",
            appearance={"message_box": "warning"},
        )
        submit_widget_state()
        exit()
    cluster_options_by_object.append(cluster_options)

if len(cluster_options_by_object) > 1:
    shared_cluster_keys = set.intersection(*[set(keys) for keys in cluster_options_by_object])
    cluster_options = [
        key
        for key in cluster_options_by_object[0]
        if key in shared_cluster_keys
    ]
else:
    cluster_options = cluster_options_by_object[0]

if not cluster_options:
    w_text_output(
        content="No shared WT/GE cluster metadata columns were found for side-by-side plotting.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

notebook_palettes = await get_notebook_palettes()
palette = w_select(
    label="Palette",
    key="count_barplot_palette",
    default=DEFAULT_CATEGORICAL_PALETTE_NAME,
    options=get_palette_selector_options(notebook_palettes),
    appearance={
        "help_text": "Use a palette saved from the H5 Viewer or fall back to the default palette."
    },
)

feature_default = choose_default_option(
    feature_options,
    preferred="ACTB",
    fallback=feature_options[0],
)
feature_select = w_select(
    label="Feature",
    key="count_feature",
    default=feature_default,
    options=tuple(feature_options),
    appearance={"help_text": "Select a feature present in the plotted object(s)."},
)

cluster_select = w_select(
    label="Cluster",
    key="umi_cluster_key",
    default=choose_group_default(cluster_options),
    options=tuple(cluster_options),
    appearance={"help_text": "Select the cluster labels used for grouping."},
)

w_row(items=[feature_select, cluster_select, palette])

feature = feature_select.value
cluster_key = cluster_select.value


def count_source_message(object_name, counts_source, adata, feature):
    if counts_source == "counts":
        return f"{object_name}: `{object_name}.layers['counts']`"
    if counts_source == "raw":
        return f"{object_name}: `{object_name}.raw.X`"

    raw_missing_feature = (
        getattr(adata, "raw", None) is not None
        and adata.raw.X is not None
        and feature not in pd.Index(adata.raw.var_names).astype(str)
    )
    if raw_missing_feature:
        return (
            f"{object_name}: `{object_name}.X` "
            f"(`{feature}` is not present in `.raw.var_names`)"
        )
    return f"{object_name}: `{object_name}.X`"


summary_frames = []
source_messages = []
for obj in count_objects:
    label = obj["label"]
    selected_adata = obj["adata"]
    selected_object_name = obj["object_name"]

    if feature not in selected_adata.var_names:
        w_text_output(
            content=f"`{feature}` was not found in `{selected_object_name}.var_names`.",
            appearance={"message_box": "warning"},
        )
        submit_widget_state()
        exit()

    if cluster_key not in selected_adata.obs:
        w_text_output(
            content=f"`{cluster_key}` was not found in `{selected_object_name}.obs`.",
            appearance={"message_box": "warning"},
        )
        submit_widget_state()
        exit()

    counts_matrix, counts_source, counts_feature_names = get_counts_matrix_for_feature(
        selected_adata,
        feature,
    )

    if feature not in counts_feature_names:
        w_text_output(
            content=f"`{feature}` was not found in `{selected_object_name}` count source.",
            appearance={"message_box": "warning"},
        )
        submit_widget_state()
        exit()

    feature_idx = list(counts_feature_names).index(feature)
    feature_counts = matrix_column_to_array(counts_matrix, feature_idx)
    cluster_labels = selected_adata.obs[cluster_key].astype(str)

    plot_df = pd.DataFrame({
        "object": label,
        "cluster": cluster_labels.values,
        "counts": feature_counts,
    })
    summary_df = (
        plot_df
        .groupby(["object", "cluster"], sort=False)
        .agg(
            mean_counts_per_spot=("counts", "mean"),
            total_counts=("counts", "sum"),
            n_spots=("counts", "size"),
            pct_spots_detected=("counts", lambda x: float((x > 0).mean() * 100)),
        )
        .reset_index()
    )
    summary_df["object"] = summary_df["object"].astype(str)
    summary_df["cluster"] = summary_df["cluster"].astype(str)
    summary_frames.append(summary_df)
    source_messages.append(
        count_source_message(
            selected_object_name,
            counts_source,
            selected_adata,
            feature,
        )
    )

summary_df = pd.concat(summary_frames, ignore_index=True)
cluster_order = sort_group_categories(summary_df["cluster"].unique().tolist())
object_order = [obj["label"] for obj in count_objects]
selected_colors = get_selected_palette_colors(
    notebook_palettes,
    palette.value,
    fallback_colors=DEFAULT_H5_CATEGORICAL_PALETTE,
)
color_map = build_discrete_color_map(cluster_order, selected_colors)

fig = px.bar(
    summary_df,
    x="mean_counts_per_spot",
    y="cluster",
    color="cluster",
    facet_col="object",
    facet_col_spacing=0.08,
    color_discrete_map=color_map,
    category_orders={
        "cluster": cluster_order,
        "object": object_order,
    },
    orientation="h",
    hover_data={
        "object": True,
        "cluster": True,
        "mean_counts_per_spot": ":.3f",
        "total_counts": ":.0f",
        "n_spots": True,
        "pct_spots_detected": ":.1f",
    },
)
fig.update_layout(
    title=f"{feature} mean counts per spot by {cluster_key}",
    xaxis_title="Mean counts per spot",
    yaxis_title=cluster_key,
    plot_bgcolor="rgba(0,0,0,0)",
    showlegend=False,
    height=max(450, 180 + 24 * len(cluster_order)),
)
fig.for_each_annotation(lambda annotation: annotation.update(text=annotation.text.split("=")[-1]))
fig.update_yaxes(categoryorder="array", categoryarray=cluster_order, matches=None)
fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor="lightgrey")

w_text_output(
    content="Count sources: " + "; ".join(source_messages) + ".",
    appearance={"message_box": "info"},
)
w_plot(source=fig)

show_table = w_checkbox(
    label="Display summary table",
    key="umi_barplot_table",
    default=False,
)

if show_table.value:
    w_table(
        label=f"{feature} count summary by {cluster_key}",
        source=summary_df,
    )
