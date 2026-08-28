# Configs

Analysis parameters (which predictors, which variant class, allele-count cutoff, colors/labels)
live here as YAML, not hardcoded in notebooks. To add a predictor, extend the relevant
`config_correlations.yaml` — the tool set is read at runtime by
[`utils/variant_filtering.py`](../utils/variant_filtering.py) (`load_config`, `load_variant_class`,
`scan_variants`, `pick_annos`), so a new entry in the config is enough; no notebook needs editing.

The three setups score different underlying column sets (Genebass single-variant betas,
All-by-All META effect sizes, individual-level UK Biobank annotations), so their configs are
**not interchangeable** even where filenames match — `configs/genebass/config_correlations.yaml`
names a column `polyphen` where `configs/all_x_all/` and `configs/ukbb/` name it `polyphen_score`,
and `configs/genebass/` carries `msa_pairformer_llr` and pLDDT/TED-domain variant classes the
others don't have. Each setup's notebooks read only their own subdirectory.

## `configs/genebass/`

| File | Used by |
|---|---|
| `config_correlations.yaml` | predictor set, colors, labels, direction — every `genebass/analysis/*` notebook |
| `config_variant_classes.yaml` | variant-class filter expressions (missense, pLoF, pLDDT strata, TED domain strata, …) |
| `config_categories.yaml` | annotation category → color mapping used by the correlation heatmaps |

## `configs/all_x_all/`

| File | Used by |
|---|---|
| `config_correlations.yaml` | predictor set for the correlation heatmap notebook |
| `config_variant_classes.yaml` | variant-class filter expressions |

## `configs/ukbb/`

| File | Used by |
|---|---|
| `config_correlations.yaml` | predictor set, colors, labels, direction, covariates |
| `config_variant_classes.yaml` | variant-class filter expressions |
| `config_odds.yaml` | odds-ratio / top-N variant applets |
| `config_odds_categories.yaml` | category → color mapping for odds-based plots |
| `config_expAssays.yaml` | experimental-assay (SGE, MaveDB) predictor set |
