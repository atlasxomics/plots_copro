w_text_output(content="## H5 Viewer: `ge_glue.h5ad`")

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
    w_table(label="Metadata (ge_glue.obs)", source=adata_ge.obs)
