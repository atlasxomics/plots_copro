w_text_output(content="## Save RNA H5 Data")

new_data_signal()
rna_save_button = w_button(label="Save RNA H5 Data")
if rna_save_button.value:
    if adata_rna is None or rna_path is None:
        w_text_output(
            content="No RNA data is loaded. Load data via the **Select Data** tab first.",
            appearance={"message_box": "warning"},
        )
        submit_widget_state()
    else:
        w_text_output(
            content=f"Writing RNA data to disk (`{rna_path.name()}`)...",
            key="rna_writing_message",
            appearance={"message_box": "info"},
        )
        submit_widget_state()
        try:
            adata_rna.write(Path(rna_path.name()))
        except Exception as e:
            w_text_output(
                content=f"Write to disk failed: {e}",
                key="rna_writing_failed",
                appearance={"message_box": "warning"},
            )
            submit_widget_state()
            exit()

        w_text_output(
            content="Uploading RNA data to Latch Data...",
            key="rna_upload_message",
            appearance={"message_box": "info"},
        )
        submit_widget_state()
        try:
            rna_path.upload_from(Path(rna_path.name()))
        except Exception as e:
            w_text_output(
                content=f"Upload failed: {e}",
                key="rna_upload_failed",
                appearance={"message_box": "warning"},
            )
            submit_widget_state()
            exit()

        w_text_output(
            content=f"Upload complete — RNA data saved to `{rna_object_name}.h5ad` in Latch Data.",
            key="rna_upload_success",
            appearance={"message_box": "success"},
        )
        submit_widget_state()
