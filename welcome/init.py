import anndata
import math
import warnings
import numpy as np
import pandas as pd
import plotly.express as px
import scanpy as sc
import scipy.sparse as sp

from anndata import AnnData
from pathlib import Path
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist
from typing import List

from lplots import palettes, submit_widget_state
from lplots.reactive import Signal
from lplots.widgets.checkbox import w_checkbox
from lplots.widgets.h5 import w_h5
from lplots.widgets.igv import w_igv, IGVOptions
from lplots.widgets.ldata import w_ldata_picker
from lplots.widgets.plot import w_plot
from lplots.widgets.row import w_row
from lplots.widgets.select import w_select
from lplots.widgets.table import w_table
from lplots.widgets.text import w_text_input, w_text_output

# Suppress benign pandas display-formatting RuntimeWarning emitted when very
# large float values are cast for rendering; it does not affect plots or data.
warnings.filterwarnings(
    "ignore",
    message="overflow encountered in cast",
    category=RuntimeWarning,
)

w_text_output(content="# **Spatial Single Cell Coprofiling Report**")
w_text_output(content="""

This notebook provides interactive viewers and starter figure generation for
`atx_glue` outputs from spatial whole transcriptome plus ATAC/CUT&Tag
coprofiling experiments.

""")

DEFAULT_H5_CATEGORICAL_PALETTE = [
    "#C33530", "#282E66", "#43884A", "#7E2F8A", "#E48341",
    "#FAE64D", "#8E9ECD", "#B570A8", "#E0C3DA", "#9FD3E2",
    "#96C56C", "#E38180", "#9584B9", "#C25434", "#63B9A8",
    "#694D99", "#33707A", "#731F1C", "#D0A970", "#3D3D3D",
]
DEFAULT_CATEGORICAL_PALETTE_NAME = "Default H5 Viewer Palette"
MAX_DEFAULT_CATEGORIES = 30
MAX_PLOT_CATEGORIES = 100
OBS_ID_KEYWORDS = (
    "barcode",
    "barcodes",
    "cell_id",
    "cellid",
    "cell_name",
    "cellname",
    "spot_id",
    "spotid",
    "uuid",
)
NA_KEYS = ["barcode", "on_off", "row", "col", "xcor", "ycor", "score"]


# Globals ------------------------------------------------------------------

if "new_data_signal" not in globals():
    new_data_signal = Signal(False)
if "refresh_ge_h5_signal" not in globals():
    refresh_ge_h5_signal = Signal(False)
if "refresh_rna_h5_signal" not in globals():
    refresh_rna_h5_signal = Signal(False)

adata_ge = None
adata_rna = None
ge_path = None
rna_path = None
ge_object_name = "atac_gs_copro"
rna_object_name = "rna_copro"
outputs_dir = None
coverages_dir = None
peak2gene_dir = None
coverage_tracks = []
coverage_track_groups = {}
available_genes = []
available_ge_features = []
DEFAULT_DATASET_OBSM_KEY = "spatial_offset"


# Functions ----------------------------------------------------------------

def empty_notebook_palettes():
    return {"categorical": [], "continuous": []}


async def get_notebook_palettes():
    try:
        palette_data = await palettes.get()
    except RuntimeError as e:
        if "Event loop is closed" in str(e):
            return empty_notebook_palettes()
        print(f"Unable to load notebook palettes: {e}")
        return empty_notebook_palettes()
    except Exception as e:
        print(f"Unable to load notebook palettes: {e}")
        return empty_notebook_palettes()

    if not isinstance(palette_data, dict):
        return empty_notebook_palettes()

    return {
        kind: [
            palette for palette in palette_data.get(kind, [])
            if palette.get("display_name") and palette.get("colors")
        ]
        for kind in ("categorical", "continuous")
    }


def get_palette_selector_options(
    palette_data,
    kind="categorical",
    fallback_name=DEFAULT_CATEGORICAL_PALETTE_NAME,
):
    palette_names = []
    seen = {fallback_name}

    for palette in palette_data.get(kind, []):
        display_name = palette["display_name"]
        if display_name not in seen:
            palette_names.append(display_name)
            seen.add(display_name)

    return tuple([fallback_name] + palette_names)


def get_selected_palette_colors(
    palette_data,
    selected_name,
    kind="categorical",
    fallback_colors=None,
    fallback_name=DEFAULT_CATEGORICAL_PALETTE_NAME,
):
    if fallback_colors is None:
        fallback_colors = DEFAULT_H5_CATEGORICAL_PALETTE

    if not selected_name or selected_name == fallback_name:
        return fallback_colors

    for palette in palette_data.get(kind, []):
        if palette["display_name"] == selected_name:
            return palette["colors"]

    return fallback_colors


def build_discrete_color_map(categories, colors):
    if not colors:
        colors = DEFAULT_H5_CATEGORICAL_PALETTE
    return {category: colors[i % len(colors)] for i, category in enumerate(categories)}


def create_proportion_dataframe(
    adata: anndata.AnnData,
    group_by: str,
    stack_by: str,
    return_type: str = "proportion",
) -> pd.DataFrame:
    if group_by not in adata.obs:
        raise KeyError(group_by)
    if stack_by not in adata.obs:
        raise KeyError(stack_by)

    count_df = pd.crosstab(
        adata.obs[group_by].astype(str),
        adata.obs[stack_by].astype(str),
    )

    if return_type == "proportion":
        result_df = count_df.div(count_df.sum(axis=1), axis=0).fillna(0)
    elif return_type == "counts":
        result_df = count_df
    else:
        raise ValueError("return_type must be 'proportion' or 'counts'.")

    long_df = result_df.reset_index().melt(
        id_vars=group_by,
        value_name="value",
        var_name=stack_by,
    )
    long_df.columns = ["group_by", "stack_by", "value"]

    return long_df


def rename_obs_keys(adata: anndata.AnnData) -> anndata.AnnData:
    key_map = {
        "Sample": "sample",
        "nFrags": "n_fragment",
        "Condition": "condition",
        "Clusters": "cluster",
    }

    keys = list(adata.obs.columns)
    for src, dest in key_map.items():
        if src in keys and dest not in keys:
            adata.obs[dest] = adata.obs[src]

    return adata


def reorder_obs_columns(adata, first_col="cluster"):
    if first_col not in adata.obs.columns:
        return
    new_order = [first_col] + [c for c in adata.obs.columns if c != first_col]
    adata.obs = adata.obs[new_order]


def sort_group_categories(values):
    num_vals = []
    all_numeric = True
    for value in values:
        try:
            num_vals.append(float(value))
        except (ValueError, TypeError):
            all_numeric = False
            break

    if all_numeric:
        sorted_pairs = sorted(zip(values, num_vals), key=lambda x: x[1])
        return [value for value, _ in sorted_pairs]

    return sorted(map(str, values))


def get_categorical_obs_keys(adata: anndata.AnnData) -> List[str]:
    return [
        key for key in adata.obs.columns
        if (
            pd.api.types.is_object_dtype(adata.obs[key])
            or pd.api.types.is_categorical_dtype(adata.obs[key])
            or pd.api.types.is_bool_dtype(adata.obs[key])
        )
    ]


def get_obs_category_summary(adata: anndata.AnnData) -> List[dict]:
    summary = []
    for key in get_categorical_obs_keys(adata):
        values = adata.obs[key].dropna()
        n_unique = int(values.nunique())
        key_lower = key.lower()
        is_id_like = any(token in key_lower for token in OBS_ID_KEYWORDS)
        plot_ok = not is_id_like and 2 <= n_unique <= MAX_PLOT_CATEGORIES

        summary.append({
            "key": key,
            "n_unique": n_unique,
            "is_id_like": is_id_like,
            "plot_ok": plot_ok,
        })

    return summary


def get_groupable_obs_keys(
    adata: anndata.AnnData,
    *,
    max_categories: int = MAX_PLOT_CATEGORIES,
) -> List[str]:
    return [
        row["key"]
        for row in get_obs_category_summary(adata)
        if row["plot_ok"] and row["n_unique"] <= max_categories
    ]


def choose_default_option(options, preferred=None, fallback=None):
    options = list(options or [])
    if preferred in options:
        return preferred
    if fallback in options:
        return fallback
    return options[0] if options else None


def choose_group_default(
    options,
    preferred=("sg_leiden_merged", "sg_leiden", "cluster", "condition", "sample"),
    fallback=None,
):
    options = list(options or [])
    for key in preferred:
        if key in options:
            return key
    if fallback in options:
        return fallback
    return options[0] if options else None


def get_cluster_keys(adata: anndata.AnnData) -> List[str]:
    preferred = [
        key for key in ["sg_leiden_merged", "sg_leiden", "cluster"]
        if key in adata.obs and adata.obs[key].nunique(dropna=True) > 1
    ]
    other = [
        key for key in get_groupable_obs_keys(adata)
        if key not in preferred
        and any(token in key.lower() for token in ["cluster", "leiden", "louvain"])
    ]
    return preferred + other


def process_matrix_layout(
    adata_all,
    n_rows: int,
    n_cols: int,
    sample_key: str = "sample",
    spatial_key: str = "spatial",
    new_obsm_key: str = "X_dataset",
    tile_spacing: float = 100.0,
    sample_order_mode: str = "original",
    condition_key: str = "condition",
):
    if sample_key not in adata_all.obs or spatial_key not in adata_all.obsm:
        return

    if sample_order_mode == "original":
        samples = list(pd.unique(adata_all.obs[sample_key]))
    elif sample_order_mode == "sample":
        samples = sorted(adata_all.obs[sample_key].astype(str).unique().tolist())
    elif sample_order_mode == "condition":
        if condition_key not in adata_all.obs:
            raise ValueError(
                f"Cannot use sample_order_mode='condition' because "
                f"`adata.obs['{condition_key}']` is missing."
            )
        obs_tmp = adata_all.obs[[sample_key, condition_key]].copy()
        cond_per_sample = (
            obs_tmp
            .assign(_i=np.arange(len(obs_tmp)))
            .sort_values("_i")
            .groupby(sample_key, sort=False)[condition_key]
            .first()
        )
        samples = (
            cond_per_sample.reset_index()
            .sort_values([condition_key, sample_key], kind="stable")
            [sample_key]
            .astype(str)
            .tolist()
        )
    else:
        raise ValueError(
            "sample_order_mode must be one of {'original','sample','condition'}"
        )

    n_samples = len(samples)
    total_positions = n_rows * n_cols
    if n_samples > total_positions:
        raise ValueError(
            f"Not enough grid positions ({n_rows}x{n_cols}={total_positions}) "
            f"for {n_samples} samples"
        )

    X_new = np.empty_like(adata_all.obsm[spatial_key], dtype=float)
    grid_bounds = {}
    sample_positions = {}

    for idx, sample_name in enumerate(samples):
        row = idx // n_cols
        col = idx % n_cols
        sample_positions[sample_name] = (row, col)

        mask = adata_all.obs[sample_key].astype(str) == str(sample_name)
        xspa = adata_all.obsm[spatial_key][mask]

        l_max = xspa.max(axis=0)
        l_min = xspa.min(axis=0)
        grid_bounds[(row, col)] = {
            "width": float(l_max[0] - l_min[0]),
            "height": float(l_max[1] - l_min[1]),
            "min_x": float(l_min[0]),
            "max_y": float(l_max[1]),
        }

    row_heights = [
        max(
            (
                grid_bounds[(r, c)]["height"]
                for c in range(n_cols)
                if (r, c) in grid_bounds
            ),
            default=0.0,
        )
        for r in range(n_rows)
    ]
    col_widths = [
        max(
            (
                grid_bounds[(r, c)]["width"]
                for r in range(n_rows)
                if (r, c) in grid_bounds
            ),
            default=0.0,
        )
        for c in range(n_cols)
    ]

    row_y_offsets = [0.0]
    for i in range(n_rows - 1):
        row_y_offsets.append(row_y_offsets[-1] - row_heights[i] - tile_spacing)

    col_x_offsets = [0.0]
    for i in range(n_cols - 1):
        col_x_offsets.append(col_x_offsets[-1] + col_widths[i] + tile_spacing)

    for sample_name in samples:
        row, col = sample_positions[sample_name]
        mask = adata_all.obs[sample_key].astype(str) == str(sample_name)
        xspa = adata_all.obsm[spatial_key][mask].copy().astype(float)

        bounds = grid_bounds[(row, col)]
        xspa[:, 0] += col_x_offsets[col] - bounds["min_x"]
        xspa[:, 1] += row_y_offsets[row] - bounds["max_y"]
        X_new[mask] = xspa

    adata_all.obsm[new_obsm_key] = X_new


def prepare_adata_for_viewer(adata: anndata.AnnData) -> anndata.AnnData:
    adata = rename_obs_keys(adata)

    for col in ["n_fragment", "n_counts", "total_counts"]:
        if col in adata.obs.columns:
            adata.obs[col] = pd.to_numeric(adata.obs[col], errors="ignore")

    for group in get_groupable_obs_keys(adata):
        if adata.obs[group].dtype != object:
            adata.obs[group] = adata.obs[group].astype(str)

    if (
        "sample" in adata.obs
        and "spatial" in adata.obsm
        and DEFAULT_DATASET_OBSM_KEY not in adata.obsm
    ):
        n_samples = adata.obs["sample"].nunique()
        if n_samples > 0:
            n_cols = min(2, max(1, n_samples))
            n_rows = math.ceil(n_samples / n_cols)
            process_matrix_layout(
                adata,
                n_rows=n_rows,
                n_cols=n_cols,
                tile_spacing=300,
                new_obsm_key=DEFAULT_DATASET_OBSM_KEY,
                sample_order_mode="sample",
            )

    if "cluster" not in adata.obs and "sg_leiden_merged" in adata.obs:
        adata.obs["cluster"] = adata.obs["sg_leiden_merged"].astype(str)
    elif "cluster" not in adata.obs and "sg_leiden" in adata.obs:
        adata.obs["cluster"] = adata.obs["sg_leiden"].astype(str)

    reorder_obs_columns(adata)
    if "orig.ident" in adata.obs.columns:
        adata.obs = adata.obs.drop(columns=["orig.ident"])

    return adata


def get_rna_counts_matrix(adata: anndata.AnnData):
    if "counts" in adata.layers:
        return adata.layers["counts"], "counts"
    if getattr(adata, "raw", None) is not None and adata.raw.X is not None:
        return adata.raw.X, "raw"
    return adata.X, "X"


def get_counts_matrix_for_feature(adata: anndata.AnnData, feature: str):
    var_names = pd.Index(adata.var_names).astype(str)
    if "counts" in adata.layers:
        return adata.layers["counts"], "counts", var_names

    if getattr(adata, "raw", None) is not None and adata.raw.X is not None:
        raw_var_names = pd.Index(adata.raw.var_names).astype(str)
        if feature in raw_var_names:
            return adata.raw.X, "raw", raw_var_names

    return adata.X, "X", var_names


def matrix_column_to_array(matrix, index: int) -> np.ndarray:
    col = matrix[:, index]
    if sp.issparse(col):
        col = col.toarray()
    return np.asarray(col).ravel()


def get_feature_values_for_plot(adata: anndata.AnnData, feature: str):
    counts_matrix, counts_source, counts_feature_names = get_counts_matrix_for_feature(
        adata,
        feature,
    )
    if feature not in counts_feature_names:
        raise KeyError(feature)

    feature_idx = list(counts_feature_names).index(feature)
    return matrix_column_to_array(counts_matrix, feature_idx), counts_source


def create_violin_data(adata: anndata.AnnData, group_by: str, plot_data: str, data_type="obs"):
    if group_by not in adata.obs:
        raise KeyError(group_by)

    if data_type == "obs":
        if plot_data not in adata.obs:
            raise KeyError(plot_data)
        values = pd.to_numeric(adata.obs[plot_data], errors="coerce")
        source = "obs"
    elif data_type == "feature":
        values, source = get_feature_values_for_plot(adata, plot_data)
    else:
        raise ValueError("data_type must be either 'obs' or 'feature'.")

    df = pd.DataFrame({
        "group": adata.obs[group_by].astype(str).values,
        "value": np.asarray(values).ravel(),
    })
    df = df[pd.notna(df["value"])]

    return df, source


def collect_coverage_tracks(root_dir):
    tracks = []
    pending = [root_dir]
    root_prefix = str(root_dir.path).rstrip("/") + "/"
    while pending:
        current = pending.pop(0)
        try:
            children = list(current.iterdir())
        except Exception:
            continue

        for child in children:
            try:
                if child.is_dir():
                    pending.append(child)
                    continue
            except Exception:
                pass

            path = str(child.path)
            lower_path = path.lower()
            track_name = path[len(root_prefix):] if path.startswith(root_prefix) else path.split("/")[-1]
            if lower_path.endswith(".bw") or lower_path.endswith(".bigwig"):
                tracks.append({
                    "name": track_name,
                    "type": "wig",
                    "format": "bigwig",
                    "url": path,
                    "autoscale": True,
                    "visibilityWindow": 1000000000,
                })
            elif lower_path.endswith(".bedgraph") or lower_path.endswith(".bedgraph.gz"):
                tracks.append({
                    "name": track_name,
                    "type": "wig",
                    "format": "bedgraph",
                    "url": path,
                    "autoscale": True,
                    "visibilityWindow": 1000000000,
                })

    return sorted(tracks, key=lambda x: x["name"])


def collect_coverage_track_groups(root_dir):
    groups = {}
    group_labels = {
        "ATAC_cluster_coverages": "atac_cluster",
        "CoPro_cluster_coverages": "copro_cluster",
        "RNA_cluster_coverages": "rna_cluster",
        "atac_cluster_coverages": "atac_cluster",
        "glue_cluster_coverages": "copro_cluster",
        "rna_cluster_coverages": "rna_cluster",
        "sample_coverages": "sample",
        "condition_coverages": "condition",
        "metadata_coverages": "metadata",
    }
    for track in collect_coverage_tracks(root_dir):
        name_parts = track["name"].split("/")
        top_folder = name_parts[0] if len(name_parts) > 1 else "coverages"
        group_name = group_labels.get(top_folder, top_folder)
        display_track = dict(track)
        display_track["name"] = "/".join(name_parts[1:]) if len(name_parts) > 1 else name_parts[-1]
        groups.setdefault(group_name, []).append(display_track)

    preferred_order = [
        "copro_cluster",
        "atac_cluster",
        "rna_cluster",
        "sample",
        "condition",
        "metadata",
    ]
    ordered_groups = {
        group: groups[group]
        for group in preferred_order
        if group in groups
    }
    for group in sorted(groups):
        if group not in ordered_groups:
            ordered_groups[group] = groups[group]

    return {
        group: sorted(tracks, key=lambda x: x["name"])
        for group, tracks in ordered_groups.items()
    }


def collect_peak2gene_tracks(root_dir):
    tracks = []
    pending = [root_dir]
    root_prefix = str(root_dir.path).rstrip("/") + "/"
    while pending:
        current = pending.pop(0)
        try:
            children = list(current.iterdir())
        except Exception:
            continue

        for child in children:
            try:
                if child.is_dir():
                    pending.append(child)
                    continue
            except Exception:
                pass

            path = str(child.path)
            if not path.lower().endswith(".bedpe"):
                continue

            track_name = path[len(root_prefix):] if path.startswith(root_prefix) else path.split("/")[-1]
            tracks.append({
                "name": track_name,
                "type": "interact",
                "format": "bedpe",
                "url": path,
                "arcType": "proportional",
                "arcOrientation": "UP",
                "color": "rgb(33, 113, 181)",
                "alpha": 0.35,
                "height": 120,
                "visibilityWindow": 10000000,
            })

    return sorted(tracks, key=lambda x: x["name"])


def collect_peak2gene_track_groups(root_dir):
    tracks = collect_peak2gene_tracks(root_dir)
    if not tracks:
        return {}

    direct_tracks = []
    gene_tracks = []
    other_tracks = []
    for track in tracks:
        display_track = dict(track)
        name_parts = track["name"].split("/")
        if name_parts[0] == "bedpe":
            display_track["name"] = "/".join(name_parts[1:]) if len(name_parts) > 1 else name_parts[-1]
            direct_tracks.append(display_track)
        elif name_parts[0] == "genes_of_interest":
            display_track["name"] = "/".join(name_parts[1:]) if len(name_parts) > 1 else name_parts[-1]
            gene_tracks.append(display_track)
        else:
            other_tracks.append(display_track)

    groups = {}
    if direct_tracks:
        groups["peak2gene"] = sorted(direct_tracks, key=lambda x: x["name"])
    if gene_tracks:
        groups["peak2gene_genes"] = sorted(gene_tracks, key=lambda x: x["name"])
    if other_tracks:
        groups["peak2gene_other"] = sorted(other_tracks, key=lambda x: x["name"])

    return groups


def add_peak2gene_overlays_to_coverage_groups(coverage_groups, peak2gene_groups):
    overlay_tracks = peak2gene_groups.get("peak2gene", [])
    if not coverage_groups or not overlay_tracks:
        return coverage_groups

    merged_groups = {}
    for group, tracks in coverage_groups.items():
        merged_tracks = list(tracks)
        for track in overlay_tracks:
            overlay_track = dict(track)
            overlay_track["name"] = f"Peak2Gene / {track['name']}"
            merged_tracks.append(overlay_track)
        merged_groups[group] = merged_tracks

    return merged_groups


def cluster_marker_to_dataframe(raw_value, key):
    """Convert an AnnData `.uns` value into a DataFrame."""
    if isinstance(raw_value, pd.DataFrame):
        return raw_value.copy()
    try:
        return pd.DataFrame(raw_value)
    except Exception as err:
        raise ValueError(f"Could not convert `adata.uns['{key}']` to a table: {err}")


def _first_existing_column(df, candidates):
    for col in candidates:
        if col in df.columns:
            return col
    return None


def normalize_ranked_genes_per_cluster(raw_value, key, groupby="cluster"):
    """Normalize ArchR-style ranked gene tables into Scanpy-like DEG columns."""
    df = cluster_marker_to_dataframe(raw_value, key)
    if df.empty:
        raise ValueError(f"`adata.uns['{key}']` is empty.")

    group_col = _first_existing_column(df, [groupby, "cluster", "group_name", "group"])
    gene_col = _first_existing_column(df, ["names", "name", "gene", "feature", "features"])
    logfc_col = _first_existing_column(
        df,
        ["logfoldchanges", "Log2FC", "log2FC", "log2fc", "avg_log2FC"],
    )
    pval_col = _first_existing_column(
        df,
        ["pvals", "Pval", "p_val", "pval", "PValue", "p_value"],
    )
    padj_col = _first_existing_column(
        df,
        ["pvals_adj", "FDR", "fdr", "p_val_adj", "pval_adj", "padj"],
    )

    missing = []
    if group_col is None:
        missing.append("cluster/group")
    if gene_col is None:
        missing.append("gene/name")
    if logfc_col is None:
        missing.append("Log2FC/logfoldchanges")
    if missing:
        raise ValueError(
            f"`adata.uns['{key}']` is missing required column(s): "
            + ", ".join(missing)
        )

    out = pd.DataFrame({
        groupby: df[group_col].astype(str),
        "names": df[gene_col].astype(str),
        "logfoldchanges": pd.to_numeric(df[logfc_col], errors="coerce"),
    })

    if pval_col is not None:
        out["pvals"] = pd.to_numeric(df[pval_col], errors="coerce")
    elif padj_col is not None:
        out["pvals"] = pd.to_numeric(df[padj_col], errors="coerce")
    else:
        out["pvals"] = 0.0

    if padj_col is not None:
        out["pvals_adj"] = pd.to_numeric(df[padj_col], errors="coerce")

    out = out.dropna(subset=[groupby, "names", "logfoldchanges", "pvals"])
    if out.empty:
        raise ValueError(f"`adata.uns['{key}']` has no usable marker rows.")

    return out


def parse_cluster_marker_int(value, default, label, minimum=1):
    try:
        parsed = int(str(value).strip())
    except Exception:
        raise ValueError(f"{label} must be an integer.")
    if parsed < minimum:
        raise ValueError(f"{label} must be at least {minimum}.")
    return parsed


def parse_cluster_marker_float(value, default, label, minimum=None):
    raw = str(value).strip()
    if raw == "":
        return default
    try:
        parsed = float(raw)
    except Exception:
        raise ValueError(f"{label} must be a number.")
    if minimum is not None and parsed < minimum:
        raise ValueError(f"{label} must be at least {minimum}.")
    return parsed


def cluster_marker_sort_key(cluster):
    cluster_str = str(cluster)
    return (0, int(cluster_str)) if cluster_str.isdigit() else (1, cluster_str)


def cluster_marker_zscore_heatmap(values_df):
    return values_df.apply(
        lambda col: (col - col.mean()) / (col.std() + 1e-9),
        axis=0,
    ).clip(-3, 3)


def get_cached_cluster_marker_heatmap(
    adata,
    heatmap_key="cluster_marker_heatmap",
):
    if heatmap_key not in adata.uns:
        return None
    heatmap_df = cluster_marker_to_dataframe(adata.uns[heatmap_key], heatmap_key)
    heatmap_df = heatmap_df.apply(pd.to_numeric, errors="coerce")
    heatmap_df = heatmap_df.dropna(axis=0, how="all").dropna(axis=1, how="all")
    heatmap_df.index = heatmap_df.index.map(str)
    heatmap_df.columns = heatmap_df.columns.map(str)
    if heatmap_df.empty:
        return None
    return heatmap_df


def choose_heatmap_layer(adata):
    for layer in ["log1p", "lognorm", "normalized", "counts"]:
        if layer in adata.layers:
            return layer
    return None


def compute_cluster_marker_heatmap_from_degs(
    adata,
    deg_df,
    top_n,
    pval_col,
    pval_cutoff,
    log2fc_cutoff,
    order_mode,
    user_order,
    deg_key="cluster_marker_degs",
    groupby="cluster",
    layer=None,
):
    if layer is None:
        layer = choose_heatmap_layer(adata)
    deg_group_col = groupby
    if deg_group_col not in deg_df.columns:
        deg_group_col = _first_existing_column(
            deg_df,
            ["cluster", "group", "group_name"],
        )

    required_cols = {"names", "logfoldchanges", pval_col}
    if deg_group_col is None:
        required_cols.add(groupby)
    missing = sorted(required_cols - set(deg_df.columns))
    if missing:
        raise ValueError(f"`adata.uns['{deg_key}']` is missing columns: {missing}.")
    if groupby not in adata.obs:
        raise ValueError(f"This plot requires `adata.obs['{groupby}']`.")
    if layer is not None and layer not in adata.layers:
        raise ValueError(f"This plot requires `adata.layers['{layer}']`.")

    deg_df = deg_df.copy()
    if deg_group_col != groupby:
        deg_df[groupby] = deg_df[deg_group_col]
    deg_df[groupby] = deg_df[groupby].astype(str)
    deg_df["names"] = deg_df["names"].astype(str)
    deg_df["logfoldchanges"] = pd.to_numeric(
        deg_df["logfoldchanges"],
        errors="coerce",
    )
    deg_df[pval_col] = pd.to_numeric(deg_df[pval_col], errors="coerce")
    deg_df = deg_df.dropna(subset=[groupby, "names", "logfoldchanges", pval_col])
    deg_df = deg_df[
        (deg_df[pval_col] <= pval_cutoff)
        & (deg_df["logfoldchanges"] >= log2fc_cutoff)
    ]

    if deg_df.empty:
        raise ValueError("No DEGs remain after the selected filters.")

    clusters = sorted(deg_df[groupby].unique().tolist(), key=cluster_marker_sort_key)
    top_genes_per_cluster = {}
    for cluster in clusters:
        cluster_df = deg_df[deg_df[groupby] == cluster].head(top_n)
        top_genes_per_cluster[cluster] = cluster_df["names"].tolist()

    seen = set()
    all_top_genes = []
    for cluster in clusters:
        for gene in top_genes_per_cluster[cluster]:
            if gene not in seen:
                all_top_genes.append(gene)
                seen.add(gene)

    if len(all_top_genes) == 0:
        raise ValueError("No marker genes are available for the selected filters.")

    gene_idx = adata.var_names.get_indexer(all_top_genes)
    valid = gene_idx >= 0
    genes = [gene for gene, keep in zip(all_top_genes, valid) if keep]
    gene_idx = gene_idx[valid]
    if len(genes) == 0:
        raise ValueError(
            "None of the selected marker genes are present in `adata.var_names`."
        )

    obs_clusters = adata.obs[groupby].astype(str)
    X = adata.layers[layer] if layer is not None else adata.X
    mean_expr = pd.DataFrame(index=clusters, columns=genes, dtype=float)
    for cluster in clusters:
        mask = obs_clusters == cluster
        if int(mask.sum()) == 0:
            mean_expr.loc[cluster] = np.nan
            continue
        sub = X[mask.to_numpy(), :][:, gene_idx]
        if sp.issparse(sub):
            sub = sub.toarray()
        mean_expr.loc[cluster] = np.asarray(sub).mean(axis=0)

    mean_expr = mean_expr.dropna(axis=0, how="all").dropna(axis=1, how="all")
    if mean_expr.empty:
        raise ValueError("Unable to compute mean expression for the selected genes.")

    clusters = mean_expr.index.tolist()
    if order_mode == "DEG similarity" and len(clusters) > 1:
        values = mean_expr.to_numpy(dtype=float)
        scaled = (values - values.mean(axis=0)) / (values.std(axis=0) + 1e-9)
        cluster_order = [
            clusters[i]
            for i in leaves_list(linkage(pdist(scaled), method="ward"))
        ]
    elif order_mode == "user-selected":
        requested = [
            item.strip()
            for item in str(user_order).replace(";", ",").split(",")
            if item.strip()
        ]
        requested = [cluster for cluster in requested if cluster in clusters]
        cluster_order = requested + [
            cluster for cluster in clusters if cluster not in requested
        ]
    else:
        cluster_order = sorted(clusters, key=cluster_marker_sort_key)

    valid_gene_set = set(mean_expr.columns)
    seen_ordered = set()
    ordered_genes = []
    for cluster in cluster_order:
        for gene in top_genes_per_cluster.get(cluster, []):
            if gene in valid_gene_set and gene not in seen_ordered:
                ordered_genes.append(gene)
                seen_ordered.add(gene)

    if len(ordered_genes) == 0:
        raise ValueError("No ordered marker genes are available for the heatmap.")

    heatmap_values = mean_expr.loc[cluster_order, ordered_genes].astype(float)
    return cluster_marker_zscore_heatmap(heatmap_values)
