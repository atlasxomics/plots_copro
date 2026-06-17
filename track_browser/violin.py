
w_text_output(content="""
## Violin Plot

Plot numeric spot metadata or feature counts across groups from the selected WT or GE object.
""")

object_options = []
object_map = {}
feature_map = {}
object_file_map = {}
if adata_rna is not None:
    object_options.append("WT")
    object_map["WT"] = adata_rna
    feature_map["WT"] = available_genes
    object_file_map["WT"] = rna_object_name
if adata_ge is not None:
    object_options.append("GE")
    object_map["GE"] = adata_ge
    feature_map["GE"] = available_ge_features
    object_file_map["GE"] = ge_object_name

object_select = w_select(
    label="Object",
    key="violin_object",
    default="WT" if "WT" in object_options else object_options[0],
    options=tuple(object_options),
    appearance={"help_text": "Choose the WT/RNA or GE AnnData object."},
)

selected_label = object_select.value
selected_adata = object_map[selected_label]
selected_features = list(feature_map[selected_label])
selected_object_name = object_file_map[selected_label]

group_options = get_groupable_obs_keys(selected_adata)
if not group_options:
    w_text_output(
        content=f"No groupable metadata columns were found in `{selected_object_name}.h5ad`.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

numeric_metadata = [
    key for key in selected_adata.obs_keys()
    if key not in NA_KEYS and pd.api.types.is_numeric_dtype(selected_adata.obs[key])
]
value_options = numeric_metadata + selected_features
if not value_options:
    w_text_output(
        content=f"No numeric metadata or features were found in `{selected_object_name}.h5ad`.",
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

value_default = choose_default_option(
    value_options,
    preferred="n_counts",
    fallback=value_options[0],
)
value_select = w_select(
    label="Value",
    key=f"violin_value_{selected_label.lower()}",
    default=value_default,
    options=tuple(value_options),
    appearance={"help_text": "Select numeric metadata or a feature to plot."},
)

group_select = w_select(
    label="Group",
    key=f"violin_group_{selected_label.lower()}",
    default=choose_group_default(group_options),
    options=tuple(group_options),
    appearance={"help_text": "Select the grouping shown on the x-axis."},
)

plot_type = w_select(
    label="Plot type",
    key="violin_plot_type",
    default="box",
    options=("box", "violin"),
    appearance={"help_text": "Use box for faster rendering on large datasets."},
)

w_row(items=[object_select, value_select, group_select, plot_type, palette])

plot_value = value_select.value
group_key = group_select.value
data_type = "obs" if plot_value in numeric_metadata else "feature"

try:
    violin_df, value_source = create_violin_data(
        selected_adata,
        group_key,
        plot_value,
        data_type=data_type,
    )
except KeyError:
    w_text_output(
        content=f"`{plot_value}` or `{group_key}` was not found in `{selected_object_name}.h5ad`.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

if violin_df.empty:
    w_text_output(
        content=f"No non-missing values were found for `{plot_value}`.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

violin_categories = sort_group_categories(violin_df["group"].unique().tolist())
selected_colors = get_selected_palette_colors(
    notebook_palettes,
    palette.value,
    fallback_colors=DEFAULT_H5_CATEGORICAL_PALETTE,
)
violin_color_map = build_discrete_color_map(violin_categories, selected_colors)

if plot_type.value == "box":
    violin_fig = px.box(
        violin_df,
        x="group",
        y="value",
        points=False,
        color="group",
        color_discrete_map=violin_color_map,
        category_orders={"group": violin_categories},
    )
elif plot_type.value == "violin":
    violin_fig = px.violin(
        violin_df,
        x="group",
        y="value",
        box=True,
        points=False,
        color="group",
        color_discrete_map=violin_color_map,
        category_orders={"group": violin_categories},
    )
else:
    w_text_output(
        content="Plot type not recognized.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

violin_fig.update_layout(
    title=f"{selected_label} {plot_value} distribution by {group_key}",
    xaxis_title=group_key,
    yaxis_title=plot_value,
    plot_bgcolor="rgba(0,0,0,0)",
    showlegend=False,
)
violin_fig.update_xaxes(
    showgrid=False,
    categoryorder="array",
    categoryarray=violin_categories,
)
violin_fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor="lightgrey")

if value_source == "obs":
    source_message = f"Value source: `{selected_object_name}.obs['{plot_value}']`."
elif value_source == "counts":
    source_message = f"Value source: `{selected_object_name}.layers['counts']`."
elif value_source == "raw":
    source_message = f"Value source: `{selected_object_name}.raw.X`."
else:
    source_message = f"Value source: `{selected_object_name}.X`."

w_text_output(content=source_message, appearance={"message_box": "info"})
w_plot(source=violin_fig)

show_table = w_checkbox(
    label="Display violin data",
    key="violin_data_table",
    default=False,
)

if show_table.value:
    w_table(
        label=f"{selected_label} {plot_value} by {group_key}",
        source=violin_df,
    )
