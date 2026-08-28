# Configs

Analysis parameters (which predictors, which variant class, allele-count cutoff, colors/labels)
live here as YAML, not hardcoded in notebooks. To add a predictor, extend the relevant
`config_correlations.yaml` — the tool set is read at runtime by
[`utils/variant_filtering.py`](../utils/variant_filtering.py) (`load_config`, `load_variant_class`,
`scan_variants`, `pick_annos`), so a new entry in the config is enough; no notebook needs editing.
`pick_annos` intersects the config's predictor names against whatever columns actually exist in
a given master table, so an entry that doesn't apply to a particular dataset is silently skipped
rather than erroring — the mechanism that makes the shared configs below safe to share.

## Root — shared by `genebass/` and `all_x_all/`

`genebass` and `all_x_all` score the same predictors against the same variant-class definitions
(both reconstruct the same benchmark, just from different summary-statistics sources), so their
configs are **one shared set** at the root, sized to the union of what their notebooks actually
use — not a merge of every config either setup ever had, most of which belonged to notebooks that
didn't make it into this repo.

| File | Used by |
|---|---|
| `config_correlations.yaml` | predictor set, colors, labels, direction — every `genebass/analysis/*` notebook and both `all_x_all/` notebooks |
| `config_variant_classes.yaml` | variant-class filter expressions (missense, pLoF, pLDDT strata, TED domain strata, …) |
| `config_categories.yaml` | annotation category → color mapping; genebass-only consumer (`mean_phenotype_master_table.ipynb`), but colocated here since there's no naming clash and no reason to split it out |

One real naming difference survives the merge: genebass's annotation parquet names two columns
`polyphen`/`sift`, while all_x_all's names the same predictors `polyphen_score`/`sift_score`.
`config_correlations.yaml` lists **both** spellings under `missense` — each setup's master table
is only known to have one of the two, so `pick_annos`'s schema intersection picks the right one
automatically and drops the other, rather than either setup silently losing the predictor to a
name mismatch. This assumes genebass's `ANNO_PATH` never carries `polyphen_score`/`sift_score`
alongside `polyphen`/`sift` — unverified against the actual upstream schema. If it ever does,
both spellings would survive `KEEP_COLS` into the master table and `pick_annos` would return
both, doubling those two rows under an identical "PolyPhen2"/"SIFT" label in any notebook that
uses them — visually obvious if it happens, but worth a one-time schema check before relying on
`genebass/utils/02_create_master_table.ipynb`'s next rebuild.

Everything else that differed between the two setups' old per-setup configs (all_x_all's
`indel_lengths`, `protein_domains`, `encode_*`, `clinvar` categories, and genebass's
`missense_interface*`/`enhancer_encode`/`tf_encode` variant classes) belonged to notebooks that
were **not** carried into this repo — dropped rather than merged in, so the shared file stays
sized to what's actually here instead of accumulating dead categories.

## `configs/ukbb/`

The primary UK Biobank RAP pipeline scores a different, larger set of annotation columns (its
own `covariates` block, several predictor categories the summary-statistics reconstructions have
no equivalent for) and isn't a strict superset or subset of the root configs, so it keeps its own
directory rather than being folded in.

| File | Used by |
|---|---|
| `config_correlations.yaml` | predictor set, colors, labels, direction, covariates |
| `config_variant_classes.yaml` | variant-class filter expressions |
| `config_odds.yaml` | odds-ratio / top-N variant applets |
| `config_odds_categories.yaml` | category → color mapping for odds-based plots |
| `config_expAssays.yaml` | experimental-assay (SGE, MaveDB) predictor set |
