new_data_signal()

w_text_output(content="""

# Cluster Marker Heatmap

Create marker-feature heatmaps from ranked gene tables or SpatialGlue cluster
DEG results stored in the loaded GE and WT/RNA AnnData objects.

""")

data_options = []
data_map = {}
if adata_ge is not None and isinstance(adata_ge, AnnData):
    data_options.append("GE")
    data_map["GE"] = {
        "adata": adata_ge,
        "deg_key": "ranked_genes_per_cluster",
        "deg_params_key": None,
        "heatmap_key": None,
        "heatmap_params_key": None,
        "mode": "ranked",
        "preferred_cluster_key": "cluster",
    }
    if (
        "stagate_cluster_marker_degs" in adata_ge.uns
        or "stagate_cluster_marker_heatmap" in adata_ge.uns
    ):
        data_options.append("GE-SpatialGlue")
        data_map["GE-SpatialGlue"] = {
            "adata": adata_ge,
            "deg_key": "stagate_cluster_marker_degs",
            "deg_params_key": "stagate_cluster_marker_degs_params",
            "heatmap_key": "stagate_cluster_marker_heatmap",
            "heatmap_params_key": "stagate_cluster_marker_heatmap_params",
            "mode": "deg",
            "preferred_cluster_key": "sg_leiden_merged",
        }
if adata_rna is not None and isinstance(adata_rna, AnnData):
    data_options.append("WT/RNA")
    data_map["WT/RNA"] = {
        "adata": adata_rna,
        "deg_key": "cluster_marker_degs",
        "deg_params_key": "cluster_marker_degs_params",
        "heatmap_key": "cluster_marker_heatmap",
        "heatmap_params_key": "cluster_marker_heatmap_params",
        "mode": "deg",
        "preferred_cluster_key": None,
    }
    if (
        "stagate_cluster_marker_degs" in adata_rna.uns
        or "stagate_cluster_marker_heatmap" in adata_rna.uns
    ):
        data_options.append("RNA-SpatialGlue")
        data_map["RNA-SpatialGlue"] = {
            "adata": adata_rna,
            "deg_key": "stagate_cluster_marker_degs",
            "deg_params_key": "stagate_cluster_marker_degs_params",
            "heatmap_key": "stagate_cluster_marker_heatmap",
            "heatmap_params_key": "stagate_cluster_marker_heatmap_params",
            "mode": "deg",
            "preferred_cluster_key": "sg_leiden_merged",
        }

if not data_options:
    w_text_output(
        content="No data loaded...",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

heatmap_data_select = w_select(
    key="cluster_marker_heatmap_data_source",
    label="Data",
    default=(
        "RNA-SpatialGlue"
        if "RNA-SpatialGlue" in data_options
        else ("WT/RNA" if "WT/RNA" in data_options else data_options[0])
    ),
    options=tuple(data_options),
    appearance={"help_text": "Choose which loaded AnnData object to plot."},
)

data_config = data_map[heatmap_data_select.value]
hm_adata = data_config["adata"]
data_label = heatmap_data_select.value

deg_key = data_config["deg_key"]
deg_params_key = data_config["deg_params_key"]
heatmap_key = data_config["heatmap_key"]
heatmap_params_key = data_config["heatmap_params_key"]

notebook_palettes = await get_notebook_palettes()

hm_palette = w_select(
    key="cluster_marker_heatmap_palette",
    label="Colorscale",
    default="Default Cluster Marker Heatmap Colorscale",
    options=get_palette_selector_options(
        notebook_palettes,
        kind="continuous",
        fallback_name="Default Cluster Marker Heatmap Colorscale",
    ),
    appearance={
        "help_text": "Use a continuous palette saved from the H5 Viewer or fall back to the default marker heatmap colors."
    },
)

heatmap_params = hm_adata.uns.get(heatmap_params_key, {}) if heatmap_params_key else {}
deg_params = hm_adata.uns.get(deg_params_key, {}) if deg_params_key else {}
default_top_n = 50
if isinstance(heatmap_params, dict) and heatmap_params.get("marker_top_n") is not None:
    default_top_n = int(heatmap_params.get("marker_top_n"))
elif isinstance(deg_params, dict) and deg_params.get("marker_top_n") is not None:
    default_top_n = int(deg_params.get("marker_top_n"))

top_n_input = w_text_input(
    key="cluster_marker_heatmap_top_n",
    label="Top features per cluster",
    default=str(default_top_n),
    appearance={"help_text": "Number of top DEG rows per cluster to include."},
)

preferred_cluster_key = data_config["preferred_cluster_key"]
stored_cluster_key = None
if isinstance(deg_params, dict):
    stored_cluster_key = deg_params.get("groupby") or deg_params.get("cluster_key")

if stored_cluster_key is not None and stored_cluster_key in hm_adata.obs:
    cluster_key = stored_cluster_key
elif preferred_cluster_key is not None and preferred_cluster_key in hm_adata.obs:
    cluster_key = preferred_cluster_key
else:
    cluster_key = choose_group_default(get_cluster_keys(hm_adata))

if deg_key in hm_adata.uns:
    if data_config["mode"] == "ranked":
        if cluster_key is None:
            cluster_key = "cluster"
        candidate_degs = normalize_ranked_genes_per_cluster(
            hm_adata.uns[deg_key],
            deg_key,
            groupby=cluster_key,
        )
    else:
        candidate_degs = cluster_marker_to_dataframe(hm_adata.uns[deg_key], deg_key)
else:
    candidate_degs = None

pval_options = [
    col for col in ["pvals", "pvals_adj"]
    if candidate_degs is not None and col in candidate_degs.columns
]
if not pval_options:
    pval_options = ["pvals"]

pval_col_select = w_select(
    key="cluster_marker_heatmap_pval_col",
    label="P-value column",
    default=pval_options[0],
    options=tuple(pval_options),
)

pval_cutoff_input = w_text_input(
    key="cluster_marker_heatmap_pval_cutoff",
    label="P-value cutoff",
    default="0.05",
)

log2fc_cutoff_input = w_text_input(
    key="cluster_marker_heatmap_log2fc_cutoff",
    label="Log2FC cutoff",
    default="0.25",
)

order_select = w_select(
    key="cluster_marker_heatmap_order",
    label="Cluster order",
    default="DEG similarity",
    options=("DEG similarity", "numeric", "user-selected"),
)

default_order = ""
if cluster_key is not None:
    default_order = ", ".join(
        sorted(
            hm_adata.obs[cluster_key].astype(str).unique().tolist(),
            key=cluster_marker_sort_key,
        )
    )

user_order_input = w_text_input(
    key="cluster_marker_heatmap_user_order",
    label="User cluster order",
    default=default_order,
    appearance={
        "help_text": "Comma-separated cluster order. Used only when Cluster order is user-selected."
    },
)

w_row(items=[heatmap_data_select, hm_palette, top_n_input, pval_col_select])
w_row(items=[pval_cutoff_input, log2fc_cutoff_input])
w_row(items=[user_order_input, order_select])

try:
    if candidate_degs is not None:
        if cluster_key is None:
            raise ValueError("No cluster-like metadata column was found.")
        top_n = parse_cluster_marker_int(
            top_n_input.value,
            default_top_n,
            "Top features per cluster",
        )
        pval_cutoff = parse_cluster_marker_float(
            pval_cutoff_input.value,
            0.05,
            "P-value cutoff",
            minimum=0.0,
        )
        log2fc_cutoff = parse_cluster_marker_float(
            log2fc_cutoff_input.value,
            0.25,
            "Log2FC cutoff",
        )
        marker_heatmap_df = compute_cluster_marker_heatmap_from_degs(
            hm_adata,
            candidate_degs,
            top_n=top_n,
            pval_col=pval_col_select.value,
            pval_cutoff=pval_cutoff,
            log2fc_cutoff=log2fc_cutoff,
            order_mode=order_select.value,
            user_order=user_order_input.value,
            deg_key=deg_key,
            groupby=cluster_key,
        )
        title = f"{data_label}: Top {top_n} Marker Features per Cluster"
    else:
        marker_heatmap_df = (
            get_cached_cluster_marker_heatmap(hm_adata, heatmap_key=heatmap_key)
            if heatmap_key is not None
            else None
        )
        if marker_heatmap_df is None:
            if data_config["mode"] == "ranked":
                raise ValueError(
                    f"No ranked gene table (`adata.uns['{deg_key}']`) was found in GE."
                )
            else:
                raise ValueError(
                    f"No dynamic DEG table (`adata.uns['{deg_key}']`) or cached heatmap "
                    f"(`adata.uns['{heatmap_key}']`) was found in {data_label}."
                )
        w_text_output(
            content=(
                "Differential stats not found, defaulting to cached heatmap. "
                "Some options will not affect the plot."
            ),
            appearance={"message_box": "warning"},
        )
        title = f"{data_label}: Cached Cluster Marker Heatmap"
except Exception as e:
    w_text_output(
        content=str(e),
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

if isinstance(deg_params, dict) and deg_params:
    fixed_filters = []
    if deg_params.get("excluded_prefixes") is not None:
        fixed_filters.append(f"excluded prefixes: {deg_params.get('excluded_prefixes')}")
    if deg_params.get("expression_layer") is not None:
        fixed_filters.append(f"expression layer: {deg_params.get('expression_layer')}")
    if fixed_filters:
        w_text_output(
            content="Stored DEG table filters: " + "; ".join(fixed_filters),
            appearance={"message_box": "info"},
        )

feature_labels = marker_heatmap_df.columns.tolist()
label_every_n = max(1, math.ceil(len(feature_labels) / 80))
visible_feature_labels = [
    feature
    for i, feature in enumerate(feature_labels)
    if i % label_every_n == 0
]

cluster_marker_heatmap = px.imshow(
    marker_heatmap_df,
    color_continuous_scale=get_selected_palette_colors(
        notebook_palettes,
        hm_palette.value,
        kind="continuous",
        fallback_colors="RdYlBu_r",
        fallback_name="Default Cluster Marker Heatmap Colorscale",
    ),
    aspect="auto",
    origin="lower",
    zmin=-3,
    zmax=3,
)

cluster_marker_heatmap.update_layout(
    title=title,
    xaxis_title="Marker feature",
    yaxis_title="Cluster",
    coloraxis_colorbar=dict(
        title="Z-score",
        title_side="right",
    ),
)

cluster_marker_heatmap.update_xaxes(
    side="bottom",
    tickmode="array",
    tickvals=visible_feature_labels,
    ticktext=visible_feature_labels,
    tickangle=90,
)
cluster_marker_heatmap.update_yaxes(
    autorange="reversed",
    tickmode="array",
    tickvals=marker_heatmap_df.index.tolist(),
    ticktext=marker_heatmap_df.index.tolist(),
)

w_plot(source=cluster_marker_heatmap)

w_table(
    label="Cluster marker heatmap data",
    source=marker_heatmap_df,
)
