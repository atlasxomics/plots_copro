import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Numpy emits "overflow encountered in cast" when very large float values are
# downcast during display formatting. Silence it locally; the warnings-module
# filter in Init does not persist across cell executions, but numpy's seterr
# (thread-local) does for this cell's execution path.
np.seterr(over="ignore")

w_text_output(content="""
## Violin Plot

Compare numeric spot metadata or feature values across groups for both the RNA
and ATAC objects, stacked vertically. The selected browser Group is used for
the x-axis of each panel.
""")

new_data_signal()
if adata_rna is None or adata_ge is None:
    w_text_output(
        content="No data loaded...",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

violin_objects = []
if adata_rna is not None:
    violin_objects.append({
        "label": "RNA",
        "adata": adata_rna,
        "features": list(available_genes),
        "object_name": rna_object_name,
    })
if adata_ge is not None:
    violin_objects.append({
        "label": "ATAC",
        "adata": adata_ge,
        "features": list(available_ge_features),
        "object_name": ge_object_name,
    })

numeric_meta_by_label = {}
value_option_sets = []
for obj in violin_objects:
    a = obj["adata"]
    numeric_metadata = [
        key for key in a.obs.columns
        if key not in NA_KEYS and pd.api.types.is_numeric_dtype(a.obs[key])
    ]
    numeric_meta_by_label[obj["label"]] = set(numeric_metadata)
    value_option_sets.append(set(numeric_metadata) | set(obj["features"]))

if len(violin_objects) > 1:
    shared_values = set.intersection(*value_option_sets)
else:
    shared_values = value_option_sets[0]

# Preserve a sensible order using the first object's listing.
first_adata = violin_objects[0]["adata"]
first_numeric = [
    key for key in first_adata.obs.columns
    if key not in NA_KEYS and pd.api.types.is_numeric_dtype(first_adata.obs[key])
]
ordered_value_options = [
    v for v in (first_numeric + list(violin_objects[0]["features"]))
    if v in shared_values
]

if not ordered_value_options:
    w_text_output(
        content="No shared numeric metadata or features were found across the RNA and ATAC objects.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()


def normalize_group_key(value):
    normalized = "".join(ch for ch in str(value).lower() if ch.isalnum())
    return normalized[:-1] if normalized.endswith("s") else normalized


def first_existing_group_key(adatas, selected_group, candidates):
    for candidate in candidates:
        if all(candidate in adata.obs for adata in adatas):
            return candidate

    shared_obs = set(adatas[0].obs.columns)
    for adata in adatas[1:]:
        shared_obs &= set(adata.obs.columns)

    selected_norm = normalize_group_key(selected_group)
    for obs_key in shared_obs:
        if normalize_group_key(obs_key) == selected_norm:
            return obs_key

    return None


selected_browser_group = coverages_group.value
violin_group_candidates = {
    "copro_cluster": ["CoPro clusters", "sg_clusters", "sg_leiden_merged", "sg_leiden", "cluster"],
    "atac_cluster": ["ATAC_cluster", "cluster"],
    "rna_cluster": ["WT_cluster", "cluster"],
    "sample": ["sample"],
    "condition": ["condition"],
}
group_key = first_existing_group_key(
    [obj["adata"] for obj in violin_objects],
    selected_browser_group,
    violin_group_candidates.get(selected_browser_group, [selected_browser_group]),
)

if group_key is None:
    w_text_output(
        content=(
            f"The selected browser group `{selected_browser_group}` does not map "
            "to a shared RNA/ATAC metadata column for the violin plot."
        ),
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

notebook_palettes = await get_notebook_palettes()
palette = w_select(
    label="Palette",
    key="violin_palette",
    default=DEFAULT_CATEGORICAL_PALETTE_NAME,
    options=get_palette_selector_options(notebook_palettes),
    appearance={
        "help_text": "Use a palette saved from the H5 Viewer or fall back to the default palette."
    },
)

value_select = w_select(
    label="Value",
    key="violin_value",
    default=choose_default_option(
        ordered_value_options,
        preferred="n_counts",
        fallback=ordered_value_options[0],
    ),
    options=tuple(ordered_value_options),
    appearance={"help_text": "Numeric metadata or feature plotted for both objects."},
)

plot_type = w_select(
    label="Plot type",
    key="violin_plot_type",
    default="box",
    options=("box", "violin"),
    appearance={"help_text": "Use box for faster rendering on large datasets."},
)

plot_value = value_select.value

violin_frames = []
for obj in violin_objects:
    data_type = "obs" if plot_value in numeric_meta_by_label[obj["label"]] else "feature"
    try:
        df, value_source = create_violin_data(
            obj["adata"],
            group_key,
            plot_value,
            data_type=data_type,
        )
    except KeyError:
        df, value_source = None, None
    violin_frames.append({
        "label": obj["label"],
        "object_name": obj["object_name"],
        "df": df,
        "source": value_source,
    })

frames_with_data = [f["df"] for f in violin_frames if f["df"] is not None and not f["df"].empty]
if not frames_with_data:
    w_text_output(
        content=f"No non-missing values were found for `{plot_value}` in either object.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

all_categories = sort_group_categories(
    pd.unique(pd.concat([df["group"] for df in frames_with_data])).tolist()
)
selected_colors = get_selected_palette_colors(
    notebook_palettes,
    palette.value,
    fallback_colors=DEFAULT_H5_CATEGORICAL_PALETTE,
)
violin_color_map = build_discrete_color_map(all_categories, selected_colors)

n_panels = len(violin_objects)
violin_fig = make_subplots(
    rows=n_panels,
    cols=1,
    subplot_titles=[obj["label"] for obj in violin_objects],
    vertical_spacing=0.25,
)

for row_idx, frame in enumerate(violin_frames, start=1):
    df = frame["df"]
    if df is None or df.empty:
        continue
    present_categories = [c for c in all_categories if c in set(df["group"])]
    for cat in present_categories:
        sub = df[df["group"] == cat]
        color = violin_color_map.get(cat)
        if plot_type.value == "violin":
            violin_fig.add_trace(
                go.Violin(
                    y=sub["value"],
                    x=[cat] * len(sub),
                    name=str(cat),
                    legendgroup=str(cat),
                    showlegend=False,
                    box_visible=True,
                    line_color=color,
                    fillcolor=color,
                    opacity=0.7,
                ),
                row=row_idx,
                col=1,
            )
        else:
            violin_fig.add_trace(
                go.Box(
                    y=sub["value"],
                    x=[cat] * len(sub),
                    name=str(cat),
                    legendgroup=str(cat),
                    showlegend=False,
                    marker_color=color,
                    boxpoints=False,
                ),
                row=row_idx,
                col=1,
            )
    violin_fig.update_yaxes(
        title_text=plot_value,
        row=row_idx,
        col=1,
        showgrid=True,
        gridwidth=1,
        gridcolor="lightgrey",
    )
    violin_fig.update_xaxes(
        title_text=group_key,
        row=row_idx,
        col=1,
        categoryorder="array",
        categoryarray=all_categories,
        showgrid=False,
    )

violin_fig.update_layout(
    title=f"{plot_value} distribution by {group_key} (RNA vs ATAC)",
    plot_bgcolor="rgba(0,0,0,0)",
    showlegend=False,
    height=max(400, 360 * n_panels),
)

source_lines = []
for frame in violin_frames:
    if frame["df"] is None:
        source_lines.append(
            f"{frame['label']}: `{plot_value}`/`{group_key}` not found in `{frame['object_name']}.h5ad`."
        )
        continue
    source = frame["source"]
    if source == "obs":
        loc = f"obs['{plot_value}']"
    elif source == "counts":
        loc = "layers['counts']"
    elif source == "raw":
        loc = "raw.X"
    else:
        loc = "X"
    source_lines.append(f"{frame['label']}: `{frame['object_name']}.{loc}`")

w_plot(source=violin_fig)

w_row(items=[value_select, plot_type, palette])

show_table = w_checkbox(
    label="Display violin data",
    key="violin_data_table",
    default=False,
)

if show_table.value:
    combined_frames = []
    for frame in violin_frames:
        if frame["df"] is not None and not frame["df"].empty:
            df = frame["df"].copy()
            df.insert(0, "object", frame["label"])
            combined_frames.append(df)
    if combined_frames:
        violin_combined_df = pd.concat(combined_frames, ignore_index=True)
        w_table(
            label=f"{plot_value} by {group_key} (RNA + ATAC)",
            source=violin_combined_df,
        )
