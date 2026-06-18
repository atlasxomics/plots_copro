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
and ATAC objects, stacked vertically. The same Value and Group selections are
applied to each object.
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
group_option_sets = []
for obj in violin_objects:
    a = obj["adata"]
    numeric_metadata = [
        key for key in a.obs.columns
        if key not in NA_KEYS and pd.api.types.is_numeric_dtype(a.obs[key])
    ]
    numeric_meta_by_label[obj["label"]] = set(numeric_metadata)
    value_option_sets.append(set(numeric_metadata) | set(obj["features"]))
    group_option_sets.append(set(get_groupable_obs_keys(a)))

if len(violin_objects) > 1:
    shared_values = set.intersection(*value_option_sets)
    shared_groups = set.intersection(*group_option_sets)
else:
    shared_values = value_option_sets[0]
    shared_groups = group_option_sets[0]

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
ordered_group_options = [
    g for g in get_groupable_obs_keys(first_adata)
    if g in shared_groups
]

if not ordered_value_options:
    w_text_output(
        content="No shared numeric metadata or features were found across the RNA and ATAC objects.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()
if not ordered_group_options:
    w_text_output(
        content="No shared groupable metadata columns were found across the RNA and ATAC objects.",
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

group_select = w_select(
    label="Group",
    key="violin_group",
    default=choose_group_default(ordered_group_options),
    options=tuple(ordered_group_options),
    appearance={"help_text": "Grouping shown on the x-axis of each panel."},
)

plot_type = w_select(
    label="Plot type",
    key="violin_plot_type",
    default="box",
    options=("box", "violin"),
    appearance={"help_text": "Use box for faster rendering on large datasets."},
)

w_row(items=[value_select, group_select, plot_type, palette])

plot_value = value_select.value
group_key = group_select.value

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

w_text_output(
    content="Value sources — " + "; ".join(source_lines),
    appearance={"message_box": "info"},
)

w_plot(source=violin_fig)

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