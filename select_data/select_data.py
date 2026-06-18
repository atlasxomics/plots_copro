w_text_output(content="""
## Select Data
<details>
<summary><i>Instructions</i></summary>

Select an `atx_glue` output directory from Latch Data. The directory should
contain reduced `rna_copro_sm.h5ad` and `atac_gs_copro_sm.h5ad` files, or the
full `rna_copro.h5ad` and `atac_gs_copro.h5ad` fallback files. It can also
include a `coverages/` subdirectory with BigWig tracks and a `peak2gene/`
subdirectory with BEDPE linkage tracks.

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

    ge_sm_matches = [f for f in children if f.name() == "atac_gs_copro_sm.h5ad"]
    ge_full_matches = [f for f in children if f.name() == "atac_gs_copro.h5ad"]
    rna_sm_matches = [f for f in children if f.name() == "rna_copro_sm.h5ad"]
    rna_full_matches = [f for f in children if f.name() == "rna_copro.h5ad"]

    if (
        len(ge_sm_matches) > 1
        or len(ge_full_matches) > 1
        or len(rna_sm_matches) > 1
        or len(rna_full_matches) > 1
    ):
        w_text_output(
            content=(
                "Found multiple AnnData files with the same expected name. "
                "Expected at most one reduced and one full file per modality."
            ),
            appearance={"message_box": "danger"},
        )
        submit_widget_state()
        exit()
    if not ge_sm_matches and not ge_full_matches:
        w_text_output(
            content=(
                "Could not find an ATAC AnnData file in the selected folder. "
                "Expected `atac_gs_copro_sm.h5ad` or `atac_gs_copro.h5ad`."
            ),
            appearance={"message_box": "danger"},
        )
        submit_widget_state()
        exit()
    if not rna_sm_matches and not rna_full_matches:
        w_text_output(
            content=(
                "Could not find an RNA AnnData file in the selected folder. "
                "Expected `rna_copro_sm.h5ad` or `rna_copro.h5ad`."
            ),
            appearance={"message_box": "danger"},
        )
        submit_widget_state()
        exit()

    ge_candidates = []
    if ge_sm_matches:
        ge_candidates.append(("atac_gs_copro_sm", ge_sm_matches[0]))
    if ge_full_matches:
        ge_candidates.append(("atac_gs_copro", ge_full_matches[0]))
    rna_candidates = []
    if rna_sm_matches:
        rna_candidates.append(("rna_copro_sm", rna_sm_matches[0]))
    if rna_full_matches:
        rna_candidates.append(("rna_copro", rna_full_matches[0]))
    ge_object_name = ge_candidates[0][0]
    ge_path = ge_candidates[0][1]
    rna_object_name = rna_candidates[0][0]
    rna_path = rna_candidates[0][1]

    w_text_output(
        content=(
            "Downloading and reading SpatialGlue AnnData files; this may take a few minutes... "
            f"Using `{rna_object_name}.h5ad` for RNA and `{ge_object_name}.h5ad` for ATAC."
        ),
        appearance={"message_box": "info"},
    )
    submit_widget_state()

    ge_load_errors = []
    adata_ge = None
    for candidate_name, candidate_path in ge_candidates:
        try:
            candidate_path.download(Path(candidate_path.name()), cache=True)
            adata_ge = sc.read_h5ad(Path(candidate_path.name()))
            ge_object_name = candidate_name
            ge_path = candidate_path
            break
        except Exception as e:
            ge_load_errors.append(f"`{candidate_name}.h5ad`: {e}")
            if candidate_name == "atac_gs_copro_sm" and len(ge_candidates) > 1:
                w_text_output(
                    content=(
                        f"Could not load `atac_gs_copro_sm.h5ad`; falling back to "
                        f"`atac_gs_copro.h5ad`. Reason: {e}"
                    ),
                    appearance={"message_box": "warning"},
                )
                submit_widget_state()

    if adata_ge is None:
        w_text_output(
            content="Error loading ATAC AnnData files: " + " | ".join(ge_load_errors),
            appearance={"message_box": "danger"},
        )
        submit_widget_state()
        exit()

    rna_load_errors = []
    adata_rna = None
    for candidate_name, candidate_path in rna_candidates:
        try:
            candidate_path.download(Path(candidate_path.name()), cache=True)
            adata_rna = sc.read_h5ad(Path(candidate_path.name()))
            rna_object_name = candidate_name
            rna_path = candidate_path
            break
        except Exception as e:
            rna_load_errors.append(f"`{candidate_name}.h5ad`: {e}")
            if candidate_name == "rna_copro_sm" and len(rna_candidates) > 1:
                w_text_output(
                    content=(
                        f"Could not load `rna_copro_sm.h5ad`; falling back to "
                        f"`rna_copro.h5ad`. Reason: {e}"
                    ),
                    appearance={"message_box": "warning"},
                )
                submit_widget_state()

    if adata_rna is None:
        w_text_output(
            content="Error loading RNA AnnData files: " + " | ".join(rna_load_errors),
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
        peak2gene_track_groups = collect_peak2gene_track_groups(peak2gene_dir)
        coverage_track_groups = add_peak2gene_overlays_to_coverage_groups(
            coverage_track_groups,
            peak2gene_track_groups,
        )
        coverage_track_groups.update(peak2gene_track_groups)

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
            f"ATAC `{ge_object_name}.h5ad` {adata_ge.n_obs} spots x {adata_ge.n_vars} features; "
            f"RNA `{rna_object_name}.h5ad` {adata_rna.n_obs} spots x {adata_rna.n_vars} genes."
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
    ge_object_name = "atac_gs_copro"
    rna_object_name = "rna_copro"
    coverage_tracks = []
    coverage_track_groups = {}
    available_genes = []
    available_ge_features = []
    refresh_ge_h5_signal(False)
    refresh_rna_h5_signal(False)
    new_data_signal(True)
