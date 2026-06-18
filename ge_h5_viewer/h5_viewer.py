new_data_signal()
if adata_ge is None:
    w_text_output(
        content="No data loaded...",
        appearance={"message_box": "warning"},
    )
    exit()


w_text_output(content=f"## H5 Viewer: `{ge_object_name}.h5ad`")

refresh_ge_h5_signal()
ge_default_obsm_key = "spatial_offset"
ge_viewer_presets = {}
if ge_default_obsm_key in adata_ge.obsm:
    ge_viewer_presets = {"default_obsm_key": ge_default_obsm_key}

for ge_default_obs_key in [
    "CoPro_cluster",
    "CoPro clusters",
    "copro_cluster",
    "sg_leiden_merged",
    "sg_leiden",
    "cluster",
]:
    if ge_default_obs_key in adata_ge.obs:
        ge_viewer_presets["default_color_by"] = {
            "type": "obs",
            "key": ge_default_obs_key,
        }
        break

w_h5(
    key="ge_h5_viewer_spatial_offset_copro_cluster",
    ann_data=adata_ge,
    viewer_presets=ge_viewer_presets,
)

ge_obs_button = w_checkbox(
    label="Display spot metadata table",
    key="ge_h5_obs_button",
    default=False,
)

if ge_obs_button.value:
    w_table(label=f"Metadata ({ge_object_name}.obs)", source=adata_ge.obs)
