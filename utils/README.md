# Shared utilities

Code used across more than one setup (`../genebass/`, `../all_x_all/`, `../ukbb/`).

## `variant_filtering.py`

The config-driven filter and predictor loader every analysis notebook imports:
`load_config`, `load_variant_class`, `scan_variants`, `pick_annos`, plus the shared
`gene_trait_tool_correlations` / `derived_schema` helpers the correlation and variance-ceiling
notebooks build on. Reads YAML from `../configs/<setup>/` — see [`../configs/README.md`](../configs/README.md)
for which file does what. This is why adding a predictor is a config edit, not a notebook edit:
extend the relevant `config_correlations.yaml` and every notebook that loads that config picks
it up.

## `REGENIE_burden_test/`

How the gene–trait associations that seed all three setups were selected: REGENIE rare-variant
burden testing (`regenie_step1*.sh`, `regenie_step2*.sh`, `regenie_both_steps_together.sh`,
`prep_files_for_REGENIE2.py`, `regenie_input_data.ipynb`) followed by
[`association_filtering/`](REGENIE_burden_test/association_filtering/) — p-value thresholding,
LOFTEE effect-direction extraction, and consolidation into the `associations_parquet` /
`config_correlations.yaml` inputs the analysis notebooks depend on. Runs on the RAP against
individual-level data; not locally runnable.

## `annotations/`

Regenerating variant annotations: `add_more_annotations.py` (bulk annotation from bigWig/FASTA/
tabix reference files), `annotate_clinvar_variants.ipynb`,
`add_ukbgym_annotations_to_bcf2parquet.ipynb`. Needs the `annotation` extra
(`uv sync --extra annotation` — see the root README) for `pyBigWig`, `pyfaidx`, `pysam`,
`biopython`, `duckdb`.

## `exp_assays/`

Experimental deep-mutational-scanning assay data used by the `other_benchmarks` notebooks in
`../genebass/analysis/` and `../ukbb/analysis/`: 
`proteingym_data.ipynb` + `annotate_proteingym_variants.ipynb` (ProteinGym; the
latter needs the `annotation` extra), `marsh_data.ipynb` (Livesey &amp; Marsh pan-protein DMS
compilation from Gen. Bio. 2025, per-UniProt predictor CSVs), `ldlr_roth_data.ipynb` (LDLR DMS, Roth et al. 2025,
supp data S1&ndash;S3).

`proteingym_uniprot_to_hgnc.csv` maps ProteinGym's UniProt entry-name filename tokens to HGNC
symbols; `build_proteingym_uniprot_to_hgnc.py` regenerates it from the UniProt REST API.

The `*_data.ipynb` notebooks each build the matching file under `../data/other_benchmarks/` and
end with an assertion that the result reproduces the copy already published on HuggingFace.
