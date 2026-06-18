w_text_output(content="""
## IGV Track Browser

Browse coverage and Peak2Gene linkage tracks from the selected output directory.
""")

new_data_signal()
if adata_rna is None:
    w_text_output(
        content="No data loaded...",
        appearance={"message_box": "warning"},
    )
    exit()

if not coverage_track_groups:
    w_text_output(
        content=(
            "No browser tracks were found. Expected BigWig/bedGraph files "
            "under `coverages/` or BEDPE files under `peak2gene/`."
        ),
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

coverage_group_options = tuple(
    group for group in coverage_track_groups.keys()
    if group != "peak2gene"
)
if not coverage_group_options:
    w_text_output(
        content="No track groups are available to display.",
        appearance={"message_box": "warning"},
    )
    submit_widget_state()
    exit()

default_coverage_group = (
    "copro_cluster"
    if "copro_cluster" in coverage_group_options
    else coverage_group_options[0]
)

coverages_group = w_select(
    label="Group",
    options=coverage_group_options,
    key="coverages_group",
    default=default_coverage_group,
    appearance={"help_text": "Select the browser coverage group and violin grouping."},
)

w_row(items=[coverages_genome, coverages_group])

selected_tracks = coverage_track_groups[coverages_group.value]

w_text_output(
    content=(
        f"Displaying {len(selected_tracks)} track(s) from "
        f"`{coverages_group.value}`."
    ),
    appearance={"message_box": "info"},
)

opts: IGVOptions = {
    "genome": coverages_genome.value,
    "tracks": selected_tracks,
}

w_igv(options=opts)
