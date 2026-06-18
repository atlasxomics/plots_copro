new_data_signal()

w_text_output(content="""
# Proportion Plot

Generate stacked bar plots where the x-axis is one metadata grouping and the
stacked segments are another metadata grouping. The y-axis can show proportions
or raw spot counts.
""")

if adata_rna is None and adata_ge is None:
    w_text_output(
        content="No data loaded...",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

object_options = []
object_map = {}
object_name_map = {}
if adata_rna is not None:
    object_options.append("RNA")
    object_map["RNA"] = adata_rna
    object_name_map["RNA"] = rna_object_name
if adata_ge is not None:
    object_options.append("ATAC")
    object_map["ATAC"] = adata_ge
    object_name_map["ATAC"] = ge_object_name

object_select = w_select(
    label="Object",
    key="proportion_object",
    default="RNA" if "RNA" in object_options else object_options[0],
    options=tuple(object_options),
    appearance={"help_text": "Choose the AnnData object to summarize."},
)

selected_label = object_select.value
selected_adata = object_map[selected_label]
selected_object_name = object_name_map[selected_label]

prop_groups = get_groupable_obs_keys(selected_adata)
if not prop_groups:
    w_text_output(
        content=f"No groupable metadata columns were found in `{selected_object_name}.h5ad`.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

notebook_palettes = await get_notebook_palettes()
categorical_palette = w_select(
    label="Palette",
    key="proportion_palette",
    default=DEFAULT_CATEGORICAL_PALETTE_NAME,
    options=get_palette_selector_options(notebook_palettes),
    appearance={
        "help_text": "Use a palette saved from the H5 Viewer or fall back to the default palette."
    },
)

group_by = w_select(
    label="Group by",
    key=f"proportion_group_by_{selected_label.lower()}",
    default=choose_default_option(prop_groups, preferred="sample", fallback=prop_groups[0]),
    options=tuple(prop_groups),
    appearance={"help_text": "Metadata grouping shown on the x-axis."},
)

stack_by = w_select(
    label="Stack by",
    key=f"proportion_stack_by_{selected_label.lower()}",
    default=choose_group_default(prop_groups),
    options=tuple(prop_groups),
    appearance={"help_text": "Metadata grouping used for stacked segments."},
)

return_type = w_select(
    label="Return type",
    key="proportion_return_type",
    default="proportion",
    options=("proportion", "counts"),
    appearance={"help_text": "Display spot counts or proportions per grouping."},
)

w_row(items=[object_select, group_by, stack_by, return_type, categorical_palette])

try:
    stacked_df = create_proportion_dataframe(
        selected_adata,
        group_by.value,
        stack_by.value,
        return_type=return_type.value,
    )
except Exception as e:
    w_text_output(
        content=f"Unable to build proportion table: {e}",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

if stacked_df.empty:
    w_text_output(
        content="No values were found for the selected metadata columns.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

stack_categories = sort_group_categories(stacked_df["stack_by"].unique().tolist())
group_categories = sort_group_categories(stacked_df["group_by"].unique().tolist())
selected_colors = get_selected_palette_colors(
    notebook_palettes,
    categorical_palette.value,
    fallback_colors=DEFAULT_H5_CATEGORICAL_PALETTE,
)
stack_color_map = build_discrete_color_map(stack_categories, selected_colors)

proportion_plot = px.bar(
    stacked_df,
    x="group_by",
    y="value",
    color="stack_by",
    barmode="stack",
    color_discrete_map=stack_color_map,
    category_orders={
        "group_by": group_categories,
        "stack_by": stack_categories,
    },
    title=f"{selected_label}: Distribution of {stack_by.value} by {group_by.value}",
)

proportion_plot.update_layout(
    xaxis_title=group_by.value,
    yaxis_title="Proportion" if return_type.value == "proportion" else "Count",
    plot_bgcolor="rgba(0,0,0,0)",
    showlegend=True,
    legend_title=stack_by.value,
)
proportion_plot.update_xaxes(
    showgrid=False,
    categoryorder="array",
    categoryarray=group_categories,
)
proportion_plot.update_yaxes(showgrid=True, gridwidth=1, gridcolor="lightgrey")

w_plot(source=proportion_plot)

show_table = w_checkbox(
    label="Display proportion data",
    key="proportion_data_table",
    default=False,
)

if show_table.value:
    w_table(
        label=f"{selected_label} {stack_by.value} by {group_by.value}",
        source=stacked_df,
    )
