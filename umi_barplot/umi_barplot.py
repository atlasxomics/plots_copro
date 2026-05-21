w_text_output(content="""
## Mean Counts Per Spot

Generate a barplot of mean counts per spot grouped by SpatialGlue cluster
for a selected RNA or GE feature.
""")

new_data_signal()
if adata_rna is None and adata_ge is None:
    w_text_output(
        content="No data loaded...",
        appearance={"message_box": "warning"},
    )
    exit()

data_options = []
data_map = {}
feature_map = {}
if adata_rna is not None:
    data_options.append("RNA")
    data_map["RNA"] = adata_rna
    feature_map["RNA"] = available_genes
if adata_ge is not None:
    data_options.append("GE")
    data_map["GE"] = adata_ge
    feature_map["GE"] = available_ge_features

count_source_select = w_select(
    label="Counts",
    key="count_barplot_source",
    default="RNA" if "RNA" in data_options else data_options[0],
    options=tuple(data_options),
    appearance={"help_text": "Choose RNA expression or GE accessibility counts."},
)

selected_label = count_source_select.value
selected_adata = data_map[selected_label]
selected_features = feature_map[selected_label]
selected_object_name = "rna_glue" if selected_label == "RNA" else "ge_glue"

cluster_options = get_cluster_keys(selected_adata)
if not cluster_options:
    w_text_output(
        content=f"No cluster-like metadata columns were found in `{selected_object_name}.h5ad`.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

if not selected_features:
    w_text_output(
        content=f"No features were found in `{selected_object_name}.h5ad`.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

notebook_palettes = await get_notebook_palettes()
palette = w_select(
    label="Palette",
    default=DEFAULT_CATEGORICAL_PALETTE_NAME,
    options=get_palette_selector_options(notebook_palettes),
    appearance={
        "help_text": "Use a palette saved from the H5 Viewer or fall back to the default palette."
    },
)

feature_default = choose_default_option(
    selected_features,
    preferred="ACTB",
    fallback=selected_features[0],
)
feature_select = w_select(
    label="Feature",
    key=f"count_feature_{selected_label.lower()}",
    default=feature_default,
    options=tuple(selected_features),
    appearance={"help_text": f"Select a feature from `{selected_object_name}.var_names`."},
)

cluster_select = w_select(
    label="Cluster",
    key="umi_cluster_key",
    default=choose_group_default(cluster_options),
    options=tuple(cluster_options),
    appearance={"help_text": "Select the cluster labels used for grouping."},
)

w_row(items=[count_source_select, feature_select, cluster_select, palette])

feature = feature_select.value
cluster_key = cluster_select.value

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

counts_matrix, counts_source = get_rna_counts_matrix(selected_adata)
if counts_source == "raw":
    counts_feature_names = pd.Index(selected_adata.raw.var_names).astype(str)
else:
    counts_feature_names = pd.Index(selected_adata.var_names).astype(str)

if feature not in counts_feature_names:
    w_text_output(
        content=f"`{feature}` was not found in the selected count source.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

feature_idx = list(counts_feature_names).index(feature)
feature_counts = matrix_column_to_array(counts_matrix, feature_idx)
cluster_labels = selected_adata.obs[cluster_key].astype(str)

plot_df = pd.DataFrame({
    "cluster": cluster_labels.values,
    "counts": feature_counts,
})
summary_df = (
    plot_df
    .groupby("cluster", sort=False)
    .agg(
        mean_counts_per_spot=("counts", "mean"),
        total_counts=("counts", "sum"),
        n_spots=("counts", "size"),
        pct_spots_detected=("counts", lambda x: float((x > 0).mean() * 100)),
    )
    .reset_index()
)
summary_df["cluster"] = summary_df["cluster"].astype(str)

cluster_order = sort_group_categories(summary_df["cluster"].tolist())
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
    color_discrete_map=color_map,
    category_orders={"cluster": cluster_order},
    orientation="h",
    hover_data={
        "cluster": True,
        "mean_counts_per_spot": ":.3f",
        "total_counts": ":.0f",
        "n_spots": True,
        "pct_spots_detected": ":.1f",
    },
)
fig.update_layout(
    title=f"{selected_label} {feature} mean counts per spot by {cluster_key}",
    xaxis_title="Mean counts per spot",
    yaxis_title=cluster_key,
    plot_bgcolor="rgba(0,0,0,0)",
    showlegend=False,
)
fig.update_yaxes(categoryorder="array", categoryarray=cluster_order)
fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor="lightgrey")

if counts_source == "counts":
    source_message = f"Count source: `{selected_object_name}.layers['counts']`."
elif counts_source == "raw":
    source_message = f"Count source: `{selected_object_name}.raw.X`."
else:
    source_message = (
        f"Count source: `{selected_object_name}.X` "
        "(no `counts` layer or `.raw` matrix was found)."
    )

w_text_output(content=source_message, appearance={"message_box": "info"})
w_plot(source=fig)

show_table = w_checkbox(
    label="Display summary table",
    key="umi_barplot_table",
    default=False,
)

if show_table.value:
    w_table(
        label=f"{selected_label} {feature} count summary by {cluster_key}",
        source=summary_df,
    )
