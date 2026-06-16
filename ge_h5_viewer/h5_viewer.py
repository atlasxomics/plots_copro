w_text_output(content=f"## H5 Viewer: `{ge_object_name}.h5ad`")

new_data_signal()
if adata_ge is None:
    w_text_output(
        content="No data loaded...",
        appearance={"message_box": "warning"},
    )
    exit()

refresh_ge_h5_signal()
w_h5(ann_data=adata_ge)

ge_obs_button = w_checkbox(
    label="Display spot metadata table",
    key="ge_h5_obs_button",
    default=False,
)

if ge_obs_button.value:
    w_table(label=f"Metadata ({ge_object_name}.obs)", source=adata_ge.obs)
