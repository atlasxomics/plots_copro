from lplots.widgets.multiselect import w_multi_select

w_text_output(content="""
## ATAC Top Marker Heatmap

Z-scored heatmap of the top 5 marker features per cluster, to highlight which
features are worth coloring the viewer by below. Toggle between the SpatialGlue
**CoPro_cluster** and the native **ATAC_cluster**, and set the cluster display
order. To inspect specific features, select them in **Genes to display**.
""")

new_data_signal()
if adata_ge is None or not isinstance(adata_ge, AnnData):
    w_text_output(content="No data loaded...", appearance={"message_box": "warning"})
    exit()

refresh_ge_h5_signal()

# cluster key -> (mode, uns key). "deg" uses ranked DEG tables; "matrix" uses a
# precomputed gene x cluster z-score matrix.
ge_hm_source_map = {
    "CoPro_cluster": ("deg", "stagate_cluster_marker_degs"),
    "ATAC_cluster": ("matrix", "genes_per_cluster_hm"),
}
ge_hm_cluster_options = tuple(k for k in ge_hm_source_map if k in adata_ge.obs)
if not ge_hm_cluster_options:
    w_text_output(
        content="No CoPro_cluster / ATAC_cluster metadata was found.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

ge_hm_top_n = 5


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


ge_hm_cluster = w_select(
    label="Cluster",
    key="ge_marker_hm_cluster",
    default="ATAC_cluster" if "ATAC_cluster" in ge_hm_cluster_options else ge_hm_cluster_options[0],
    options=ge_hm_cluster_options,
    appearance={"help_text": "Cluster labels used to group and rank marker features."},
)
ge_hm_order = w_select(
    label="Cluster order",
    key="ge_marker_hm_order",
    default="numeric",
    options=("numeric", "DEG similarity", "user-selected"),
    appearance={"help_text": "Ordering of clusters along the heatmap rows."},
)

ge_hm_cluster_key = ge_hm_cluster.value

ge_hm_row_items = [ge_hm_cluster, ge_hm_order]
ge_hm_user_order_value = ""
if ge_hm_order.value == "user-selected":
    ge_hm_default_order = ", ".join(
        sorted(adata_ge.obs[ge_hm_cluster_key].astype(str).unique().tolist(), key=cluster_marker_sort_key)
    )
    ge_hm_user_order = w_text_input(
        label="User cluster order",
        key=f"ge_marker_hm_user_order_{ge_hm_cluster_key.lower()}",
        default=ge_hm_default_order,
        appearance={"help_text": "Comma-separated order of clusters along the heatmap rows."},
    )
    ge_hm_user_order_value = ge_hm_user_order.value
    ge_hm_row_items.append(ge_hm_user_order)

ge_hm_gene_options = sorted(dict.fromkeys(adata_ge.var_names.astype(str)))
ge_hm_genes_input = w_multi_select(
    label="Genes to display",
    key="ge_marker_hm_genes_sel",
    options=ge_hm_gene_options,
    default=[],
    appearance={"help_text": "Select features to display instead of the automatic top markers. Leave empty for top markers."},
)
ge_hm_row_items.append(ge_hm_genes_input)

w_row(items=ge_hm_row_items)

ge_hm_custom_genes = [str(g) for g in (ge_hm_genes_input.value or [])]

ge_hm_mode, ge_hm_uns_key = ge_hm_source_map[ge_hm_cluster_key]
ge_hm_missing = []

try:
    if ge_hm_custom_genes:
        ge_hm_df, ge_hm_missing = _gene_list_zscore_heatmap(
            adata_ge, ge_hm_custom_genes, ge_hm_cluster_key, ge_hm_order.value, ge_hm_user_order_value
        )
    else:
        if ge_hm_uns_key not in adata_ge.uns:
            raise ValueError(
                f"No precomputed marker table (`adata.uns['{ge_hm_uns_key}']`) is available "
                f"for `{ge_hm_cluster_key}`."
            )

        if ge_hm_mode == "deg":
            ge_hm_deg_df = cluster_marker_to_dataframe(adata_ge.uns[ge_hm_uns_key], ge_hm_uns_key)
            ge_hm_adata = adata_ge
            ge_hm_df = compute_cluster_marker_heatmap_from_degs(
                ge_hm_adata,
                ge_hm_deg_df,
                top_n=ge_hm_top_n,
                pval_col="pvals" if "pvals" in ge_hm_deg_df.columns else "pvals_adj",
                pval_cutoff=0.05,
                log2fc_cutoff=0.25,
                order_mode=ge_hm_order.value,
                user_order=ge_hm_user_order_value,
                deg_key=ge_hm_uns_key,
                groupby=ge_hm_cluster_key,
            )
        else:
            ge_hm_matrix = cluster_marker_to_dataframe(adata_ge.uns[ge_hm_uns_key], ge_hm_uns_key)
            ge_hm_matrix.columns = [str(c) for c in ge_hm_matrix.columns]
            ge_hm_matrix = ge_hm_matrix.apply(pd.to_numeric, errors="coerce").dropna(axis=0, how="all")
            ge_hm_clusters = list(ge_hm_matrix.columns)
            if ge_hm_order.value == "user-selected":
                ge_hm_requested = [
                    x.strip() for x in str(ge_hm_user_order_value).replace(";", ",").split(",") if x.strip()
                ]
                ge_hm_requested = [c for c in ge_hm_requested if c in ge_hm_clusters]
                ge_hm_cluster_order = ge_hm_requested + [c for c in ge_hm_clusters if c not in ge_hm_requested]
            elif ge_hm_order.value == "DEG similarity" and len(ge_hm_clusters) > 1:
                ge_hm_cluster_order = [
                    ge_hm_clusters[i]
                    for i in leaves_list(linkage(pdist(ge_hm_matrix.T.to_numpy(dtype=float)), method="ward"))
                ]
            else:
                ge_hm_cluster_order = sorted(ge_hm_clusters, key=cluster_marker_sort_key)

            ge_hm_seen, ge_hm_genes = set(), []
            for c in ge_hm_cluster_order:
                for gsym in ge_hm_matrix[c].sort_values(ascending=False).head(ge_hm_top_n).index.tolist():
                    if gsym not in ge_hm_seen:
                        ge_hm_genes.append(gsym)
                        ge_hm_seen.add(gsym)
            ge_hm_df = ge_hm_matrix.loc[ge_hm_genes, ge_hm_cluster_order].T
            ge_hm_df.index = ge_hm_df.index.map(str)
            ge_hm_df = ge_hm_df.clip(-3, 3)
except Exception as e:
    w_text_output(content=str(e), appearance={"message_box": "warning"})
    submit_widget_state()
    exit()

if ge_hm_missing:
    w_text_output(
        content="Not found in this dataset (ignored): " + ", ".join(ge_hm_missing),
        appearance={"message_box": "warning"},
    )

ge_hm_title = (
    f"ATAC: {len(ge_hm_df.columns)} Selected Features per {ge_hm_cluster_key}"
    if ge_hm_custom_genes
    else f"ATAC: Top {ge_hm_top_n} Marker Features per {ge_hm_cluster_key}"
)

ge_hm_fig = px.imshow(
    ge_hm_df,
    color_continuous_scale="RdYlBu_r",
    aspect="auto",
    origin="lower",
    zmin=-3,
    zmax=3,
)
ge_hm_fig.update_layout(
    title=ge_hm_title,
    xaxis_title="Marker feature",
    yaxis_title=ge_hm_cluster_key,
    coloraxis_colorbar=dict(title="Z-score", title_side="right"),
    height=max(320, 90 + 30 * len(ge_hm_df.index)),
)
ge_hm_fig.update_xaxes(side="bottom", tickangle=90)
ge_hm_fig.update_yaxes(
    autorange="reversed",
    tickmode="array",
    tickvals=ge_hm_df.index.tolist(),
    ticktext=ge_hm_df.index.tolist(),
)

w_plot(source=ge_hm_fig)

ge_hm_show_table = w_checkbox(
    label="Display marker heatmap data",
    key="ge_marker_hm_table",
    default=False,
)
if ge_hm_show_table.value:
    w_table(label=f"ATAC top markers by {ge_hm_cluster_key}", source=ge_hm_df)
