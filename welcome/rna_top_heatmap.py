from lplots.widgets.multiselect import w_multi_select

w_text_output(content="""
## RNA Top Marker Heatmap

Z-scored heatmap of the top 5 marker features per cluster, to highlight which
features are worth coloring the viewer by below. Toggle between the SpatialGlue
**CoPro_cluster** and the native **RNA_cluster**, and set the cluster display
order. To inspect specific features, select them in **Genes to display**.
""")

new_data_signal()
if adata_rna is None or not isinstance(adata_rna, AnnData):
    w_text_output(content="No data loaded...", appearance={"message_box": "warning"})
    exit()

refresh_rna_h5_signal()

# cluster key -> (mode, uns key). "deg" uses ranked DEG tables. CoPro_cluster
# uses the SpatialGlue marker table; RNA_cluster uses the RNA-native
# `cluster_marker_degs` table.
rna_hm_source_map = {
    "CoPro_cluster": ("deg", "stagate_cluster_marker_degs"),
    "RNA_cluster": ("deg", "cluster_marker_degs"),
}
rna_hm_cluster_options = tuple(k for k in rna_hm_source_map if k in adata_rna.obs)
if not rna_hm_cluster_options:
    w_text_output(
        content="No CoPro_cluster / RNA_cluster metadata was found.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

rna_hm_top_n = 5


def _gene_list_zscore_heatmap(adata, genes, groupby, order_mode, user_order):
    import numpy as np
    import scipy.sparse as sp

    a = adata
    var_set = set(a.var_names)
    seen, present, missing = set(), [], []
    for g in genes:
        if g in seen:
            continue
        seen.add(g)
        (present if g in var_set else missing).append(g)
    if not present:
        raise ValueError("None of the requested genes were found: " + ", ".join(genes))
    X = read_matrix_columns(a.X, feature_column_indices(a, present))
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    expr = pd.DataFrame(X, columns=present, index=a.obs_names.astype(str))
    grp = a.obs[groupby].astype(str).to_numpy()
    mean_expr = expr.groupby(grp).mean()
    z = (mean_expr - mean_expr.mean(axis=0)) / mean_expr.std(axis=0).replace(0, np.nan)
    z = z.fillna(0.0).clip(-3, 3)
    clusters = list(z.index)
    if order_mode == "user-selected":
        req = [x.strip() for x in str(user_order).replace(";", ",").split(",") if x.strip()]
        req = [c for c in req if c in clusters]
        order = req + [c for c in clusters if c not in req]
    elif order_mode == "DEG similarity" and len(clusters) > 1:
        order = [clusters[i] for i in leaves_list(linkage(pdist(z.to_numpy(dtype=float)), method="ward"))]
    else:
        order = sorted(clusters, key=cluster_marker_sort_key)
    z = z.loc[order, present]
    z.index = z.index.map(str)
    return z, missing


rna_hm_cluster = w_select(
    label="Cluster",
    key="rna_marker_hm_cluster",
    default="RNA_cluster" if "RNA_cluster" in rna_hm_cluster_options else rna_hm_cluster_options[0],
    options=rna_hm_cluster_options,
    appearance={"help_text": "Cluster labels used to group and rank marker features."},
)
rna_hm_order = w_select(
    label="Cluster order",
    key="rna_marker_hm_order",
    default="numeric",
    options=("numeric", "DEG similarity", "user-selected"),
    appearance={"help_text": "Ordering of clusters along the heatmap rows."},
)

rna_hm_cluster_key = rna_hm_cluster.value

rna_hm_row_items = [rna_hm_cluster, rna_hm_order]
rna_hm_user_order_value = ""
if rna_hm_order.value == "user-selected":
    rna_hm_default_order = ", ".join(
        sorted(adata_rna.obs[rna_hm_cluster_key].astype(str).unique().tolist(), key=cluster_marker_sort_key)
    )
    rna_hm_user_order = w_text_input(
        label="User cluster order",
        key=f"rna_marker_hm_user_order_{rna_hm_cluster_key.lower()}",
        default=rna_hm_default_order,
        appearance={"help_text": "Comma-separated order of clusters along the heatmap rows."},
    )
    rna_hm_user_order_value = rna_hm_user_order.value
    rna_hm_row_items.append(rna_hm_user_order)

rna_hm_gene_options = sorted(dict.fromkeys(adata_rna.var_names.astype(str)))
rna_hm_genes_input = w_multi_select(
    label="Genes to display",
    key="rna_marker_hm_genes_sel",
    options=rna_hm_gene_options,
    default=[],
    appearance={"help_text": "Select features to display instead of the automatic top markers. Leave empty for top markers."},
)
rna_hm_row_items.append(rna_hm_genes_input)

w_row(items=rna_hm_row_items)

rna_hm_custom_genes = [str(g) for g in (rna_hm_genes_input.value or [])]

rna_hm_mode, rna_hm_uns_key = rna_hm_source_map[rna_hm_cluster_key]
rna_hm_missing = []

if not rna_hm_custom_genes and (rna_hm_uns_key is None or rna_hm_uns_key not in adata_rna.uns):
    w_text_output(
        content=(
            f"No precomputed marker table (`adata.uns['{rna_hm_uns_key}']`) is available "
            f"for `{rna_hm_cluster_key}`. Select features in **Genes to display** to view a heatmap."
        ),
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

try:
    if rna_hm_custom_genes:
        rna_hm_df, rna_hm_missing = _gene_list_zscore_heatmap(
            adata_rna, rna_hm_custom_genes, rna_hm_cluster_key, rna_hm_order.value, rna_hm_user_order_value
        )
    else:
        rna_hm_deg_df = cluster_marker_to_dataframe(adata_rna.uns[rna_hm_uns_key], rna_hm_uns_key)
        rna_hm_adata = adata_rna
        rna_hm_df = compute_cluster_marker_heatmap_from_degs(
            rna_hm_adata,
            rna_hm_deg_df,
            top_n=rna_hm_top_n,
            pval_col="pvals" if "pvals" in rna_hm_deg_df.columns else "pvals_adj",
            pval_cutoff=0.05,
            log2fc_cutoff=0.25,
            order_mode=rna_hm_order.value,
            user_order=rna_hm_user_order_value,
            deg_key=rna_hm_uns_key,
            groupby=rna_hm_cluster_key,
        )
except Exception as e:
    w_text_output(content=str(e), appearance={"message_box": "warning"})
    submit_widget_state()
    exit()

if rna_hm_missing:
    w_text_output(
        content="Not found in this dataset (ignored): " + ", ".join(rna_hm_missing),
        appearance={"message_box": "warning"},
    )

rna_hm_title = (
    f"RNA: {len(rna_hm_df.columns)} Selected Features per {rna_hm_cluster_key}"
    if rna_hm_custom_genes
    else f"RNA: Top {rna_hm_top_n} Marker Features per {rna_hm_cluster_key}"
)

rna_hm_fig = px.imshow(
    rna_hm_df,
    color_continuous_scale="RdYlBu_r",
    aspect="auto",
    origin="lower",
    zmin=-3,
    zmax=3,
)
rna_hm_fig.update_layout(
    title=rna_hm_title,
    xaxis_title="Marker feature",
    yaxis_title=rna_hm_cluster_key,
    coloraxis_colorbar=dict(title="Z-score", title_side="right"),
    height=max(320, 90 + 30 * len(rna_hm_df.index)),
)
rna_hm_fig.update_xaxes(side="bottom", tickangle=90)
rna_hm_fig.update_yaxes(
    autorange="reversed",
    tickmode="array",
    tickvals=rna_hm_df.index.tolist(),
    ticktext=rna_hm_df.index.tolist(),
)

w_plot(source=rna_hm_fig)

rna_hm_show_table = w_checkbox(
    label="Display marker heatmap data",
    key="rna_marker_hm_table",
    default=False,
)
if rna_hm_show_table.value:
    w_table(label=f"RNA top markers by {rna_hm_cluster_key}", source=rna_hm_df)
