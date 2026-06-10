w_text_output(content="""
## Select Data
<details>
<summary><i>Instructions</i></summary>

Select an `atx_glue` output directory from Latch Data. The directory should
contain `ge_glue.h5ad`, `rna_glue.h5ad`, and optionally a `coverages/`
subdirectory with BigWig tracks and a `peak2gene/` subdirectory with BEDPE
linkage tracks.

</details>
""")

data_path = w_ldata_picker(
    label="atx_glue output folder",
    key="data_path",
    appearance={"placeholder": "Select a glue_outs project folder"},
)

if data_path.value is not None:
    if not data_path.value.is_dir():
        w_text_output(
            content="Selected resource must be an `atx_glue` output directory.",
            appearance={"message_box": "danger"},
        )
        submit_widget_state()
        exit()

    outputs_dir = data_path.value
    children = list(outputs_dir.iterdir())

    ge_matches = [f for f in children if f.name() == "ge_glue.h5ad"]
    rna_matches = [f for f in children if f.name() == "rna_glue.h5ad"]

    if len(ge_matches) != 1:
        w_text_output(
            content="Could not find exactly one `ge_glue.h5ad` file in the selected folder.",
            appearance={"message_box": "danger"},
        )
        submit_widget_state()
        exit()
    if len(rna_matches) != 1:
        w_text_output(
            content="Could not find exactly one `rna_glue.h5ad` file in the selected folder.",
            appearance={"message_box": "danger"},
        )
        submit_widget_state()
        exit()

    ge_path = ge_matches[0]
    rna_path = rna_matches[0]

    w_text_output(
        content="Downloading and reading SpatialGlue AnnData files; this may take a few minutes...",
        appearance={"message_box": "info"},
    )
    submit_widget_state()

    try:
        ge_path.download(Path(ge_path.name()), cache=True)
        rna_path.download(Path(rna_path.name()), cache=True)
    except Exception as e:
        w_text_output(
            content=f"Error downloading input files: {e}",
            appearance={"message_box": "danger"},
        )
        submit_widget_state()
        exit()

    try:
        adata_ge = sc.read_h5ad(Path(ge_path.name()))
    except Exception as e:
        w_text_output(
            content=f"Error loading `ge_glue.h5ad`: {e}",
            appearance={"message_box": "danger"},
        )
        submit_widget_state()
        exit()

    try:
        adata_rna = sc.read_h5ad(Path(rna_path.name()))
    except Exception as e:
        w_text_output(
            content=f"Error loading `rna_glue.h5ad`: {e}",
            appearance={"message_box": "danger"},
        )
        submit_widget_state()
        exit()

    adata_ge = prepare_adata_for_viewer(adata_ge)
    adata_rna = prepare_adata_for_viewer(adata_rna)

    available_ge_features = list(adata_ge.var_names)
    available_genes = list(adata_rna.var_names)

    coverages_dir = None
    peak2gene_dir = None
    coverage_matches = [
        f for f in children
        if f.name() == "coverages" and f.is_dir()
    ]
    peak2gene_matches = [
        f for f in children
        if f.name() == "peak2gene" and f.is_dir()
    ]
    if coverage_matches:
        coverages_dir = coverage_matches[0]
        coverage_track_groups = collect_coverage_track_groups(coverages_dir)
    else:
        coverage_track_groups = {}

    if peak2gene_matches:
        peak2gene_dir = peak2gene_matches[0]
        coverage_track_groups.update(collect_peak2gene_track_groups(peak2gene_dir))

    coverage_tracks = [
        track
        for tracks in coverage_track_groups.values()
        for track in tracks
    ]

    if not coverage_track_groups:
        w_text_output(
            content=(
                "No browser tracks were found. Expected BigWig/bedGraph files "
                "under `coverages/` or BEDPE files under `peak2gene/`."
            ),
            appearance={"message_box": "warning"},
        )
        submit_widget_state()

    w_text_output(
        content=(
            "Data successfully loaded: "
            f"GE {adata_ge.n_obs} spots x {adata_ge.n_vars} features; "
            f"RNA {adata_rna.n_obs} spots x {adata_rna.n_vars} genes."
        ),
        appearance={"message_box": "success"},
    )
    submit_widget_state()

    refresh_ge_h5_signal(False)
    refresh_rna_h5_signal(False)
    new_data_signal(True)
else:
    adata_ge = None
    adata_rna = None
    ge_path = None
    rna_path = None
    outputs_dir = None
    coverages_dir = None
    peak2gene_dir = None
    coverage_tracks = []
    coverage_track_groups = {}
    available_genes = []
    available_ge_features = []
    refresh_ge_h5_signal(False)
    refresh_rna_h5_signal(False)
    new_data_signal(True)
