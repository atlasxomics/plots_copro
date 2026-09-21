import plotly.graph_objects as go
from plotly.subplots import make_subplots
from lplots.widgets.multiselect import w_multi_select

new_data_signal()

w_text_output(content="""

## Cluster Marker Heatmap — ATAC vs RNA

Z-scored marker heatmaps for the **ATAC** and **RNA** objects over a shared
marker-feature set. Rows are grouped by the chosen **Cluster grouping** (CoPro /
ATAC / RNA), and the features come from the chosen **DEG source** (top N per
cluster, filtered by the p-value and log2FC cutoffs). Use **Display** to show
one modality or both side-by-side; when both are shown they share the same
features so the modalities can be compared directly. Features absent from a
modality render as grey NA tiles. Select features in **Genes to display** to
override the automatic markers.

""")

if (
    adata_ge is None or adata_rna is None
    or not isinstance(adata_ge, AnnData) or not isinstance(adata_rna, AnnData)
):
    w_text_output(
        content="Both the ATAC and RNA objects must be loaded for the side-by-side comparison.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

refresh_ge_h5_signal()
refresh_rna_h5_signal()

RDYLBU_R = px.colors.diverging.RdYlBu[::-1]

GENE_COL_CANDIDATES = ["names", "name", "gene", "feature", "features"]
LOGFC_COL_CANDIDATES = ["logfoldchanges", "Log2FC", "log2FC", "log2fc", "avg_log2FC"]
GROUP_COL_CANDIDATES = ["cluster", "group", "group_name", "CoPro_cluster", "RNA_cluster", "ATAC_cluster"]


# Row grouping: cluster labels used to group both panels.
hm_group_options = tuple(
    k for k in ("CoPro_cluster", "ATAC_cluster", "RNA_cluster")
    if k in adata_ge.obs and k in adata_rna.obs
)
if not hm_group_options:
    w_text_output(
        content="No shared CoPro_cluster / ATAC_cluster / RNA_cluster metadata was found across both objects.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

# DEG source: differential table used to pick marker features. "deg" uses a
# ranked DEG table (respects p-value / log2FC filters); "matrix" uses a gene x
# cluster z-score matrix (filters do not apply). The SpatialGlue/CoPro marker
# table exists on both objects, so it is exposed separately for RNA-derived and
# ATAC-derived markers.
hm_deg_source_map = {}
if "stagate_cluster_marker_degs" in adata_rna.uns:
    hm_deg_source_map["CoPro (SpatialGlue) — RNA markers"] = {
        "mode": "deg", "uns_key": "stagate_cluster_marker_degs", "adata": adata_rna,
    }
if "stagate_cluster_marker_degs" in adata_ge.uns:
    hm_deg_source_map["CoPro (SpatialGlue) — ATAC markers"] = {
        "mode": "deg", "uns_key": "stagate_cluster_marker_degs", "adata": adata_ge,
    }
if "cluster_marker_degs" in adata_rna.uns:
    hm_deg_source_map["RNA"] = {
        "mode": "deg", "uns_key": "cluster_marker_degs", "adata": adata_rna,
    }
if "genes_per_cluster_hm" in adata_ge.uns:
    hm_deg_source_map["ATAC"] = {
        "mode": "matrix", "uns_key": "genes_per_cluster_hm", "adata": adata_ge,
    }
hm_deg_source_options = tuple(hm_deg_source_map.keys())

hm_default_top_n = 5


def _deg_top_genes(df, top_n, pval_col, pval_cutoff, log2fc_cutoff):
    group_col = _first_existing_column(df, GROUP_COL_CANDIDATES) or df.columns[0]
    gene_col = _first_existing_column(df, GENE_COL_CANDIDATES)
    logfc_col = _first_existing_column(df, LOGFC_COL_CANDIDATES)
    if gene_col is None:
        raise ValueError("The selected DEG table has no gene/name column.")
    d = df.copy()
    d[gene_col] = d[gene_col].astype(str)
    d[group_col] = d[group_col].astype(str)
    if logfc_col is not None:
        d[logfc_col] = pd.to_numeric(d[logfc_col], errors="coerce")
        d = d[d[logfc_col].fillna(-np.inf) >= log2fc_cutoff]
    if pval_col in d.columns:
        d[pval_col] = pd.to_numeric(d[pval_col], errors="coerce")
        d = d[d[pval_col].fillna(np.inf) <= pval_cutoff]
    d = d.dropna(subset=[gene_col])
    if d.empty:
        raise ValueError("No DEGs remain after the selected filters.")
    seen, genes = set(), []
    for cl in sorted(d[group_col].unique().tolist(), key=cluster_marker_sort_key):
        for g in d[d[group_col] == cl].head(top_n)[gene_col].tolist():
            if g not in seen:
                seen.add(g)
                genes.append(g)
    return genes


def _matrix_top_genes(df, top_n):
    m = df.copy()
    m.columns = [str(c) for c in m.columns]
    m = m.apply(pd.to_numeric, errors="coerce").dropna(axis=0, how="all")
    seen, genes = set(), []
    for c in sorted(m.columns, key=cluster_marker_sort_key):
        for g in m[c].sort_values(ascending=False).head(top_n).index.tolist():
            if str(g) not in seen:
                seen.add(str(g))
                genes.append(str(g))
    return genes


def _modality_zscore_aligned(adata, genes, groupby, cluster_order):
    # Cluster x gene z-scored matrix for `adata` aligned to the given gene and
    # cluster order. Genes/clusters absent from this modality stay NaN so they
    # render as grey NA tiles.
    import numpy as np
    import scipy.sparse as sp

    genes = [str(g) for g in genes]
    cluster_order = [str(c) for c in cluster_order]
    result = pd.DataFrame(index=cluster_order, columns=genes, dtype=float)
    if adata is None or not isinstance(adata, AnnData) or groupby not in adata.obs:
        return result
    a = adata
    var_set = set(a.var_names.astype(str))
    present = [g for g in genes if g in var_set]
    if not present:
        return result
    X = read_matrix_columns(a.X, feature_column_indices(a, present))
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    expr = pd.DataFrame(X, columns=present, index=a.obs_names.astype(str))
    grp = a.obs[groupby].astype(str).to_numpy()
    mean_expr = expr.groupby(grp).mean()
    z = (mean_expr - mean_expr.mean(axis=0)) / mean_expr.std(axis=0).replace(0, np.nan)
    z = z.clip(-3, 3)
    z.index = z.index.map(str)
    z = z.reindex(index=cluster_order, columns=present)
    for g in present:
        result[g] = z[g]
    return result


notebook_palettes = await get_notebook_palettes()

# Display toggle in its own row at the top. Defaults to RNA.
hm_display = w_select(
    key="cluster_marker_heatmap_display",
    label="Display",
    default="RNA",
    options=("Both", "ATAC", "RNA"),
    appearance={"help_text": "Show a single modality or both panels side-by-side."},
)
w_row(items=[hm_display])

hm_group = w_select(
    key="cluster_marker_heatmap_group",
    label="Cluster grouping",
    default="RNA_cluster" if "RNA_cluster" in hm_group_options else hm_group_options[0],
    options=hm_group_options,
    appearance={"help_text": "Cluster labels used to group the rows of both panels."},
)
hm_group_key = hm_group.value

# Default DEG source follows the grouping when a matching table exists.
_deg_default_by_group = {
    "CoPro_cluster": [
        "CoPro (SpatialGlue) — RNA markers",
        "CoPro (SpatialGlue) — ATAC markers",
    ],
    "RNA_cluster": ["RNA"],
    "ATAC_cluster": ["ATAC"],
}
hm_deg_default = next(
    (
        name
        for name in _deg_default_by_group.get(hm_group_key, [])
        if name in hm_deg_source_options
    ),
    hm_deg_source_options[0] if hm_deg_source_options else None,
)

hm_deg_source = w_select(
    key="cluster_marker_heatmap_deg_source",
    label="DEG source",
    default=hm_deg_default,
    options=hm_deg_source_options,
    appearance={"help_text": "Differential table used to select which marker features to display."},
)
hm_order = w_select(
    key="cluster_marker_heatmap_order",
    label="Cluster order",
    default="numeric",
    options=("numeric", "DEG similarity", "user-selected"),
    appearance={"help_text": "Ordering of clusters along the heatmap rows."},
)
hm_palette = w_select(
    key="cluster_marker_heatmap_palette",
    label="Colorscale",
    default="Default Cluster Marker Heatmap Colorscale",
    options=get_palette_selector_options(
        notebook_palettes,
        kind="continuous",
        fallback_name="Default Cluster Marker Heatmap Colorscale",
    ),
    appearance={"help_text": "Use a continuous palette saved from the H5 Viewer or the default marker colors."},
)
hm_top_n_input = w_text_input(
    key="cluster_marker_heatmap_top_n",
    label="Top features per cluster",
    default=str(hm_default_top_n),
    appearance={"help_text": "Number of top marker features per cluster to include."},
)

hm_cfg = hm_deg_source_map[hm_deg_source.value] if hm_deg_source.value in hm_deg_source_map else None
hm_is_deg_mode = hm_cfg is not None and hm_cfg["mode"] == "deg"

# Load the selected DEG table up front so the p-value column options reflect it.
hm_deg_df = None
if hm_cfg is not None and hm_cfg["uns_key"] in hm_cfg["adata"].uns:
    try:
        hm_deg_df = cluster_marker_to_dataframe(hm_cfg["adata"].uns[hm_cfg["uns_key"]], hm_cfg["uns_key"])
    except Exception:
        hm_deg_df = None

hm_pval_options = [
    col for col in ("pvals", "pvals_adj")
    if hm_is_deg_mode and hm_deg_df is not None and col in hm_deg_df.columns
] or ["pvals"]

hm_row_top = [hm_group, hm_deg_source, hm_order, hm_palette]
w_row(items=hm_row_top)

hm_pval_col_input = None
hm_pval_cutoff_input = None
hm_log2fc_cutoff_input = None
if hm_is_deg_mode:
    hm_pval_col_input = w_select(
        key="cluster_marker_heatmap_pval_col",
        label="P-value column",
        default=hm_pval_options[0],
        options=tuple(hm_pval_options),
        appearance={"help_text": "Column used for the p-value cutoff."},
    )
    hm_pval_cutoff_input = w_text_input(
        key="cluster_marker_heatmap_pval_cutoff",
        label="P-value cutoff",
        default="0.05",
        appearance={"help_text": "Keep DEGs with p-value at or below this."},
    )
    hm_log2fc_cutoff_input = w_text_input(
        key="cluster_marker_heatmap_log2fc_cutoff",
        label="Log2FC cutoff",
        default="0.25",
        appearance={"help_text": "Keep DEGs with log2 fold-change at or above this."},
    )
    w_row(items=[hm_top_n_input, hm_pval_col_input, hm_pval_cutoff_input, hm_log2fc_cutoff_input])
else:
    w_row(items=[hm_top_n_input])

hm_user_order_value = ""
hm_bottom_items = []
if hm_order.value == "user-selected":
    hm_all_clusters = sorted(
        set(adata_ge.obs[hm_group_key].astype(str).unique().tolist())
        | set(adata_rna.obs[hm_group_key].astype(str).unique().tolist()),
        key=cluster_marker_sort_key,
    )
    hm_user_order = w_text_input(
        key=f"cluster_marker_heatmap_user_order_{hm_group_key.lower()}",
        label="User cluster order",
        default=", ".join(hm_all_clusters),
        appearance={"help_text": "Comma-separated cluster order along the heatmap rows."},
    )
    hm_user_order_value = hm_user_order.value
    hm_bottom_items.append(hm_user_order)

hm_gene_options = sorted(
    dict.fromkeys(
        list(adata_ge.var_names.astype(str)) + list(adata_rna.var_names.astype(str))
    )
)
hm_genes_input = w_multi_select(
    label="Genes to display",
    key="cluster_marker_heatmap_genes_sel",
    options=hm_gene_options,
    default=[],
    appearance={"help_text": "Select features to display instead of the automatic top markers. Leave empty for top markers."},
)
hm_bottom_items.append(hm_genes_input)
w_row(items=hm_bottom_items)

hm_custom_genes = [str(g) for g in (hm_genes_input.value or [])]

# Which modality panels to show.
hm_selected_panels = ["ATAC", "RNA"] if hm_display.value == "Both" else [hm_display.value]

# 1) Select the feature set.
try:
    hm_top_n = parse_cluster_marker_int(
        hm_top_n_input.value, hm_default_top_n, "Top features per cluster"
    )
    if hm_custom_genes:
        hm_feature_order = list(dict.fromkeys(hm_custom_genes))
    elif hm_cfg is None:
        raise ValueError("No DEG source is available. Select features in Genes to display instead.")
    elif hm_is_deg_mode:
        if hm_deg_df is None:
            raise ValueError(
                f"No usable DEG table (`adata.uns['{hm_cfg['uns_key']}']`) for the selected DEG source."
            )
        hm_pval_cutoff = parse_cluster_marker_float(
            hm_pval_cutoff_input.value, 0.05, "P-value cutoff", minimum=0.0
        )
        hm_log2fc_cutoff = parse_cluster_marker_float(
            hm_log2fc_cutoff_input.value, 0.25, "Log2FC cutoff"
        )
        hm_feature_order = _deg_top_genes(
            hm_deg_df, hm_top_n, hm_pval_col_input.value, hm_pval_cutoff, hm_log2fc_cutoff
        )
    else:
        if hm_cfg["uns_key"] not in hm_cfg["adata"].uns:
            raise ValueError(
                f"No marker matrix (`adata.uns['{hm_cfg['uns_key']}']`) for the selected DEG source."
            )
        hm_matrix_df = cluster_marker_to_dataframe(
            hm_cfg["adata"].uns[hm_cfg["uns_key"]], hm_cfg["uns_key"]
        )
        hm_feature_order = _matrix_top_genes(hm_matrix_df, hm_top_n)
except Exception as e:
    w_text_output(content=str(e), appearance={"message_box": "warning"})
    submit_widget_state()
    exit()

if not hm_feature_order:
    w_text_output(
        content="No features are available to display for the current selection.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

# 2) Determine the cluster (row) order from the grouping.
hm_clusters = set(adata_ge.obs[hm_group_key].astype(str).unique().tolist()) | set(
    adata_rna.obs[hm_group_key].astype(str).unique().tolist()
)
if hm_order.value == "user-selected":
    req = [x.strip() for x in str(hm_user_order_value).replace(";", ",").split(",") if x.strip()]
    req = [c for c in req if c in hm_clusters]
    hm_cluster_order = req + [c for c in sorted(hm_clusters, key=cluster_marker_sort_key) if c not in req]
else:
    hm_cluster_order = sorted(hm_clusters, key=cluster_marker_sort_key)

# 3) Z-score each shown modality over the shared features (grey NA where missing).
hm_adata_by_panel = {"ATAC": adata_ge, "RNA": adata_rna}
hm_panel_data = {
    name: _modality_zscore_aligned(hm_adata_by_panel[name], hm_feature_order, hm_group_key, hm_cluster_order)
    for name in hm_selected_panels
}

# Optional DEG-similarity reordering of the rows using the shown panels.
if hm_order.value == "DEG similarity" and len(hm_cluster_order) > 2:
    hm_combined = pd.concat([hm_panel_data[name] for name in hm_selected_panels], axis=1).to_numpy(dtype=float)
    hm_combined = np.nan_to_num(hm_combined, nan=0.0)
    try:
        hm_leaves = leaves_list(linkage(pdist(hm_combined), method="ward"))
        hm_cluster_order = [hm_cluster_order[i] for i in hm_leaves]
        hm_panel_data = {name: df.reindex(index=hm_cluster_order) for name, df in hm_panel_data.items()}
    except Exception:
        pass

hm_panels = [(name, hm_panel_data[name]) for name in hm_selected_panels]

hm_feature_label = (
    f"{len(hm_feature_order)} selected features"
    if hm_custom_genes
    else f"top {hm_top_n} {hm_deg_source.value} DEGs"
)
hm_title_prefix = " vs ".join(hm_selected_panels) if len(hm_selected_panels) > 1 else hm_selected_panels[0]
hm_colorscale = get_selected_palette_colors(
    notebook_palettes,
    hm_palette.value,
    kind="continuous",
    fallback_colors=RDYLBU_R,
    fallback_name="Default Cluster Marker Heatmap Colorscale",
)

hm_fig = make_subplots(
    rows=1,
    cols=len(hm_panels),
    shared_yaxes=True,
    horizontal_spacing=0.12,
    subplot_titles=[name for name, _ in hm_panels],
)
for col_idx, (name, panel_df) in enumerate(hm_panels, start=1):
    hm_fig.add_trace(
        go.Heatmap(
            z=panel_df.to_numpy(dtype=float),
            x=list(panel_df.columns),
            y=list(panel_df.index),
            zmin=-3,
            zmax=3,
            colorscale=hm_colorscale,
            showscale=(col_idx == 1),
            colorbar=dict(title="Z-score", title_side="right") if col_idx == 1 else None,
            hoverongaps=False,
        ),
        row=1,
        col=col_idx,
    )
    hm_fig.update_xaxes(title_text="Marker feature", tickangle=90, row=1, col=col_idx)
    hm_fig.update_yaxes(
        title_text=hm_group_key if col_idx == 1 else None,
        autorange="reversed",
        row=1,
        col=col_idx,
    )

hm_fig.update_layout(
    title=f"{hm_title_prefix} — {hm_feature_label} per {hm_group_key}",
    plot_bgcolor="lightgrey",
    height=max(360, 140 + 30 * len(hm_cluster_order)),
)

w_plot(source=hm_fig)

hm_show_table = w_checkbox(
    label="Display heatmap data",
    key="cluster_marker_heatmap_table",
    default=False,
)
if hm_show_table.value:
    for name, panel_df in hm_panels:
        w_table(label=f"{name} ({hm_group_key})", source=panel_df)
