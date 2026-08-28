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
tabix reference files), `annotate_clinvar_variants.ipynb`, `annotate_proteingym_variants.ipynb`,
`add_ukbgym_annotations_to_bcf2parquet.ipynb`. Needs the `annotation` extra
(`uv sync --extra annotation` — see the root README) for `pyBigWig`, `pyfaidx`, `pysam`,
`biopython`, `duckdb`.

## `exp_assays/`

Experimental deep-mutational-scanning assay data used by the `other_benchmarks` notebooks in
`../genebass/analysis/` and `../ukbb/analysis/`: `extract_sge_scores.ipynb` (Saturation Genome
Editing), `mavedb_data.ipynb` (MaveDB), `proteingym_data.ipynb` (ProteinGym).
