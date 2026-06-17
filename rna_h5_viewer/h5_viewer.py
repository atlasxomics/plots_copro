w_text_output(content=f"## H5 Viewer: `{rna_object_name}.h5ad`")

new_data_signal()
if adata_rna is None:
    w_text_output(
        content="No data loaded...",
        appearance={"message_box": "warning"},
    )
    exit()

refresh_rna_h5_signal()
w_h5(ann_data=adata_rna)

rna_obs_button = w_checkbox(
    label="Display spot metadata table",
    key="rna_h5_obs_button",
    default=False,
)

if rna_obs_button.value:
    w_table(label=f"Metadata ({rna_object_name}.obs)", source=adata_rna.obs)
