new_data_signal()
if adata_rna is None:
    w_text_output(
        content="No data loaded...",
        appearance={"message_box": "warning"},
    )
    exit()

w_text_output(content=f"## H5 Viewer: `{rna_object_name}.h5ad`")

refresh_rna_h5_signal()
rna_default_obsm_key = "spatial_offset"
rna_viewer_presets = {}
if rna_default_obsm_key in adata_rna.obsm:
    rna_viewer_presets = {"default_obsm_key": rna_default_obsm_key}

for rna_default_obs_key in [
    "CoPro_cluster",
    "CoPro clusters",
    "copro_cluster",
    "sg_leiden_merged",
    "sg_leiden",
    "cluster",
]:
    if rna_default_obs_key in adata_rna.obs:
        rna_viewer_presets["default_color_by"] = {
            "type": "obs",
            "key": rna_default_obs_key,
        }
        break

w_h5(
    key="rna_h5_viewer_spatial_offset_copro_cluster",
    ann_data=adata_rna,
    viewer_presets=rna_viewer_presets,
)

rna_obs_button = w_checkbox(
    label="Display spot metadata table",
    key="rna_h5_obs_button",
    default=False,
)

if rna_obs_button.value:
    w_table(label=f"Metadata ({rna_object_name}.obs)", source=adata_rna.obs)
