w_text_output(content="""
## IGV Track Browser

Browse all coverage tracks found under the selected `coverages/` directory.
""")

new_data_signal()
if adata_rna is None:
    w_text_output(
        content="No data loaded...",
        appearance={"message_box": "warning"},
    )
    exit()

if coverages_dir is None:
    w_text_output(
        content="No `coverages/` directory was found in the selected outputs.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

if not coverage_track_groups:
    w_text_output(
        content="No BigWig or bedGraph files were found under `coverages/`.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

coverages_genome = w_select(
    label="Genome",
    options=("hg38", "mm10", "rn6"),
    key="coverages_genome",
    default="hg38",
    appearance={"help_text": "Select reference genome."},
)

coverage_group_options = tuple(coverage_track_groups.keys())
default_coverage_group = (
    "glue_cluster"
    if "glue_cluster" in coverage_track_groups
    else coverage_group_options[0]
)

coverages_group = w_select(
    label="Coverage tracks",
    options=coverage_group_options,
    key="coverages_group",
    default=default_coverage_group,
    appearance={"help_text": "Select ATAC, SpatialGlue, or RNA cluster coverage tracks."},
)

w_row(items=[coverages_genome, coverages_group])

selected_tracks = coverage_track_groups[coverages_group.value]

w_text_output(
    content=(
        f"Displaying {len(selected_tracks)} coverage track(s) from "
        f"`{coverages_group.value}`."
    ),
    appearance={"message_box": "info"},
)

opts: IGVOptions = {
    "genome": coverages_genome.value,
    "tracks": selected_tracks,
}

w_igv(options=opts)
