w_text_output(content="## Save ATAC H5 Data")

new_data_signal()
ge_save_button = w_button(label="Save ATAC H5 Data")
if ge_save_button.value:
    if adata_ge is None or ge_path is None:
        w_text_output(
            content="No ATAC data is loaded. Load data via the **Select Data** tab first.",
            appearance={"message_box": "warning"},
        )
        submit_widget_state()
    else:
        w_text_output(
            content=f"Writing ATAC data to disk (`{ge_path.name()}`)...",
            key="ge_writing_message",
            appearance={"message_box": "info"},
        )
        submit_widget_state()
        try:
            adata_ge.write(Path(ge_path.name()))
        except Exception as e:
            w_text_output(
                content=f"Write to disk failed: {e}",
                key="ge_writing_failed",
                appearance={"message_box": "warning"},
            )
            submit_widget_state()
            exit()

        w_text_output(
            content="Uploading ATAC data to Latch Data...",
            key="ge_upload_message",
            appearance={"message_box": "info"},
        )
        submit_widget_state()
        try:
            ge_path.upload_from(Path(ge_path.name()))
        except Exception as e:
            w_text_output(
                content=f"Upload failed: {e}",
                key="ge_upload_failed",
                appearance={"message_box": "warning"},
            )
            submit_widget_state()
            exit()

        w_text_output(
            content=f"Upload complete — ATAC data saved to `{ge_object_name}.h5ad` in Latch Data.",
            key="ge_upload_success",
            appearance={"message_box": "success"},
        )
        submit_widget_state()
