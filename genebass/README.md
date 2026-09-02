# Genebass Summary Statistics for UKBBGym

## Goal

Reconstruct the UKBBGym coding-variant benchmark from [Genebass](https://genebass.org/) single-variant summary statistics (394,841 UK Biobank exomes), so it can be run without individual-level data on the UK Biobank RAP. This is the setup that **runs locally**.

UKBBGym scores a variant by the **average phenotype among its carriers** (APPV) and evaluates variant effect predictors by the rank correlation between their scores and that quantity, per gene–trait association. Genebass publishes a SAIGE single-variant **BETA** instead. For a rare variant the two estimate the same thing up to a deterministic $(1-p)$ factor that cancels under any rank-based metric, so the BETA is used as the APPV proxy throughout.

Gene–trait pairs are *selected* by pLoF burden tests, but the benchmark scores variants of **every** consequence class, so all variants in the target genes are extracted, not just pLoF.

## Run this locally

You don't need the ~1 TB Genebass MatrixTable to reproduce the figures — the built master table
(`genebass_annotated.parquet`) is published in the
[`gagneurlab/ukbbgym`](https://huggingface.co/datasets/gagneurlab/ukbbgym) Hugging Face dataset,
and every notebook fetches its own inputs from there automatically on first run. Steps 1 and 2
below are only needed if you want to rebuild the master table from scratch.

1. From the repo root: `uv sync` (see the root [README](../README.md) for environment details).
2. Log in to Hugging Face with an account that has access to `gagneurlab/ukbbgym` (it's private):
   `.venv/bin/huggingface-cli login`, or export `HF_TOKEN`. No manual download or file placement
   needed — the first notebook you run pulls whatever it needs into `data/` (gitignored) itself.
   [`correlations.ipynb`](analysis/correlations.ipynb),
   [`mean_phenotype.ipynb`](analysis/mean_phenotype.ipynb) and
   [`protein_domains_correlations.ipynb`](analysis/protein_domains_correlations.ipynb)
   need only the master table; every other analysis notebook fetches an additional file or two —
   see the table below and the *needs* column in *Analysis notebooks*.
3. Open any notebook under [`analysis/`](analysis/) with the `.venv` kernel and run top to bottom,
   or run all of them at once: `genebass/run_all.sh`. It executes every notebook under
   `analysis/` in place (via `nbconvert`, against the repo's `.venv`), auto-downloading each
   notebook's inputs as it goes, and reports which notebooks failed.
4. Figures are written to `paper_figures/` at the repo root.

### Changing parameters without editing a notebook

Every notebook's parameter cell reads its defaults through `env_override()`
([`utils/variant_filtering.py`](../utils/variant_filtering.py)) — the current default when its
`UKBBGYM_<NAME>` environment variable is unset, that variable's value otherwise. `run_all.sh`
exposes the common ones as flags, applied to every notebook that has that parameter (one without
it just ignores the flag):

```bash
genebass/run_all.sh --variant-class indel --mac 10 --min-variants 50 \
                     --selected-categories "missense,conservation"
```

| Flag | Env var | Notebook default |
|---|---|---|
| `--variant-class` | `UKBBGYM_VARIANT_CLASS` | `missense` (most notebooks; `indel` for `noise_ceiling_one_region.ipynb`) |
| `--selected-categories` | `UKBBGYM_SELECTED_CATEGORIES` (comma-separated) | per-notebook — `config_correlations.yaml` categories. **Not** `mean_phenotype.ipynb` or `noise_ceiling/*` — see below, their category list means something different |
| `--mac` | `UKBBGYM_MAC` | `20` |
| `--min-variants` | `UKBBGYM_MIN_VARIANTS` | per-notebook (typically `50`–`100`) |
| `--config-file` | `UKBBGYM_CONFIG_FILE` | per-notebook, always a file under [`../configs/`](../configs/). **Not** `mean_phenotype.ipynb` — see below |
| `--master-path` | `UKBBGYM_MASTER_PATH` | `data/genebass_annotated.parquet`, auto-fetched from Hugging Face. **Not** `expAssays_all_genes_correlations_pheno.ipynb` — see below |
| `--fig-dir` | `UKBBGYM_FIG_DIR` | `paper_figures/` at the repo root |
| `--only-snps` / `--no-only-snps` | `UKBBGYM_ONLY_SNPS` | per-notebook |
| `--new-score PATH` | (sets `UKBBGYM_MASTER_PATH`/`UKBBGYM_CONFIG_FILE` itself) | none — run with the existing annotations only |
| `--new-score-category NAME` | — | `missense` |

### Scoring your own predictor: `--new-score`

To see how a new predictor compares against the existing tools in every plot, without
touching any notebook: give `run_all.sh` a parquet or CSV(.gz) with one or more score columns
plus a way to join each score to a variant — either:

- **genomic** — `chrom`, `pos` (or `position`), `ref`, `alt` (tried first; if the file also
  has `uniprot_id`, that's folded into the key too, so a score given per (variant, isoform)
  collapses onto one row per variant instead of colliding), or
- **protein** — `uniprot_id`, `position` (or `protein_position`), `aa_ref`, `aa_alt`, used if
  no genomic key is found (e.g. a per-protein score with no genomic coordinates at all, like
  an MSA-based LLR table).

A row that repeats a join key is collapsed to its first occurrence, so the merge can't
silently duplicate master-table rows. Pass `--on` (`merge_new_score.py`) to force a different
key if your file uses neither shape.

```bash
genebass/run_all.sh --new-score /path/to/my_predictor_scores.parquet
```

Before running any notebook, this:

1. Left-joins your file's new columns onto the master table by the detected key
   (`merge_new_score()` in [`utils/variant_filtering.py`](../utils/variant_filtering.py)),
   writing `data/<master file>_plus_<your file's name>.parquet`.
2. Registers each new column as a predictor under the `missense` category (override with
   `--new-score-category`) in a copy of `config_correlations.yaml`
   (`register_score_in_config()`, same file) — written to `configs/_custom_score.yaml`
   (gitignored, safe to overwrite on the next run).
3. Points every notebook at that merged table and config for the rest of the run, the same
   way `--master-path`/`--config-file` would.

`missense` is the default category because it's in every notebook's default
`selected_categories` — the new score shows up without also passing
`--selected-categories`. It only reaches notebooks that read `MASTER_PATH`/`CONFIG_FILE` off
the shared flags — not `mean_phenotype.ipynb`,
`expAssays_all_genes_correlations_pheno.ipynb`, or the `noise_ceiling/*` unassociated table,
which use the special-cased env vars below and their own master tables.

If none is specified, `run_all.sh` runs with the existing annotations exactly as before.

**Three notebooks give a same-shaped parameter its own env var instead**, because its default
isn't the same thing the flag above controls — a shared `--selected-categories`/`--config-file`/
`--master-path` would otherwise silently repoint them at the wrong config or table with no
error. Set these by exporting the env var directly; there's no flag for them:

| Notebook | Env var | What it is |
|---|---|---|
| `mean_phenotype.ipynb` | `UKBBGYM_MEAN_PHENO_CATEGORIES` | categories from `config_categories.yaml` (default `['protein_domains']`), not `config_correlations.yaml` |
| `mean_phenotype.ipynb` | `UKBBGYM_MEAN_PHENO_CONFIG_FILE` | defaults to `config_categories.yaml`, not `config_correlations.yaml` |
| `noise_ceiling/noise_ceiling_one_region.ipynb`, `noise_ceiling/noise_ceiling_TED_contrast.ipynb` | `UKBBGYM_NOISE_CEILING_CATEGORIES` | defaults to `None` (use the variant class's own `tool_categories`); the shared flag would disable that fallback |
| `other_benchmarks/expAssays_all_genes_correlations_pheno.ipynb` | `UKBBGYM_EXPASSAYS_MASTER_PATH` | defaults to `genebass_annotated_all.parquet`, the all-variants table, not the shared master table |

To run just one notebook with an override, export the env var(s) yourself and use `.venv/bin/jupyter nbconvert --execute` (or open it interactively — the env var still applies):

```bash
UKBBGYM_MAC=10 .venv/bin/jupyter nbconvert --to notebook --execute --inplace \
    genebass/analysis/correlations.ipynb
```

A few other notebook-specific knobs aren't on either list above but follow the same
`UKBBGYM_<NAME>` pattern — check each notebook's `# --- parameters` cell for the exact names
(e.g. `UKBBGYM_ONLY_CLINVAR`, `UKBBGYM_DOMAIN_TYPE`, `UKBBGYM_CI_FACTOR`, `UKBBGYM_N_BOOT`).

### `data/` inventory — fetched automatically from Hugging Face

Every `data/` input is a file in the private
[`gagneurlab/ukbbgym`](https://huggingface.co/datasets/gagneurlab/ukbbgym) dataset. You don't
place any of these by hand: each notebook's parameter cell calls `fetch_hf_data(...)`
([`utils/variant_filtering.py`](../utils/variant_filtering.py)), which checks `data/<file>`
first and, if it's missing, downloads it there via `huggingface_hub` — so the first run of a
given notebook fetches only the file(s) it actually needs, and every run after that is local.
The dataset is **private**: fetching requires a Hugging Face account with access, logged in
locally (`huggingface-cli login`, or an `HF_TOKEN`/`HUGGING_FACE_HUB_TOKEN` env var) — ask
whoever runs the lab's Hugging Face org for access if you get a 401.

| File (local `data/` path = dataset path) | Needed by |
|---|---|
| `genebass_annotated.parquet` | most `analysis/` notebooks (the master table) |
| `genebass_unassociated.parquet` | `noise_ceiling/*` (unassociated/null-trait panel, AC=1 synonymous noise estimate) |
| `genebass_annotated_all.parquet` | `other_benchmarks/expAssays_all_genes_correlations_pheno.ipynb` (all-variants variant of the master table, not restricted to the 670 benchmark genes) |
| `other_benchmarks/clinvar_annotated.parquet` | `other_benchmarks/clinvar_spearman_scatterplot.ipynb` (17,683-gene ClinVar label set; broader than the master table's `clinical_significance` column, see the notebook's own note) |
| `other_benchmarks/proteingym_snv_annotated.parquet` | `other_benchmarks/proteingym_snr.ipynb`, `other_benchmarks/proteingym_correlations.ipynb`, `other_benchmarks/expAssays_all_genes_correlations_pheno.ipynb` |
| `other_benchmarks/marsh_gen_bio_2025.parquet`, `other_benchmarks/ldlr_science_2025.parquet` | `other_benchmarks/expAssays_all_genes_correlations_pheno.ipynb` (additional experimental DMS assays) |

To point a notebook at a file you already have locally under a different name instead of
fetching, pass its path via the notebook's own `*_PATH`/`--master-path` override (see
"Changing parameters without editing a notebook" above) — `fetch_hf_data` is only the default.

## Pipeline

```
Genebass MatrixTable  ──[1]──►  genebass_betas_…parquet  ─┐
                                                          │
UKB annotation parquet  ─────────────────────────────────►├──[2]──►  master_table_…parquet  ──►  analysis/
regenie association + LOFTEE correlation files  ─────────►┘                                      (figure notebooks)
```

| Step | Notebook | Output |
|---|---|---|
| 1 | [utils/01_read_hail_sumstats.ipynb](utils/01_read_hail_sumstats.ipynb) | `genebass_betas_127phenos_allvars.parquet` — one row per (variant, phenotype) |
| 2 | [utils/02_create_master_table.ipynb](utils/02_create_master_table.ipynb) | `master_table_<date>.parquet` — one row per (variant, gene), everything joined; this is the file published on Hugging Face as `genebass_annotated.parquet` |

Step 1 needs the Genebass Hail MatrixTable (available from the Genebass authors on Google Cloud)
and a Hail/Spark environment — see *Environment setup* in the [root README](../README.md).
Step 2 needs the UKB annotation, association and LOFTEE-correlation parquets, which are only
available from the primary UK Biobank pipeline (`../ukbb/`) and are not publicly distributed.
**If you're using the published master table, skip both steps.**

Every notebook under [analysis/](analysis/) reads the master table, plus — for some notebooks —
the additional `data/` inputs listed in the inventory above. None does a further join against
the raw annotation/association files directly; anything beyond the master table is a standalone
file read wholesale (ClinVar labels, ProteinGym, the null-trait panel, …), not a join key into
`ANNO_PATH`/`ASSOC_PATH`/`CORR_PATH`.

### Path conventions

The tables below name files relative to two roots where the input data needs to be stored

| Placeholder | What it points at |
|---|---|
| `<data_root>` | the UKBBGym project data directory holding annotation, association and beta files |
| `<genebass_mt>` | the directory holding the downloaded Genebass Hail MatrixTable |
---

## Step 1 — Genebass betas

**Notebook:** [utils/01_read_hail_sumstats.ipynb](utils/01_read_hail_sumstats.ipynb)

### Input data

| File | Description |
|---|---|
| `<genebass_mt>` — `variant_results.mt` | Genebass Hail MatrixTable — 8M variants × 4,529 phenotypes, ~1 TB. Row key: `(locus, alleles)`. Entry fields: `BETA, SE, Pvalue, AF`. |
| `<data_root>/association_files/`<br>`regenie_127phenotypes_mac20_lofteeHC_EUR_correlations.parquet` | UKBBGym association file — 2,289 gene–phenotype pairs (699 genes, 121 phenotypes) to extract. `region` = Ensembl gene ID, `phenotype` = human-readable name with `_int` suffix. |
| `<data_root>/genes_info_biomart.parquet` | Ensembl ID → gene symbol mapping (Biomart export). Used to bridge the association file (Ensembl IDs) to the MT row field (`gene` = gene symbol). |
| `<data_root>/gencode/gencode.v29.annotation.gtf.gz` | GENCODE **v29** gene annotation — the version Genebass used (VEP v95). Provides each gene's `(contig, start, end)` for `locus`-interval pruning **and** the `region → symbol` map that matches the MT's `gene` field. |

### Output

`<data_root>/genebass_betas/genebass_betas_127phenos_allvars.parquet`

Flat parquet with one row per (variant, phenotype), restricted to the 2,289 association pairs. Column names match the UKBBGym pipeline:

| Column | Description |
|---|---|
| `id` | Variant ID, `chrom:pos:ref:alt` (e.g. `REDACTED_VARIANT_ID`) — joins to UKBBGym variant/score files |
| `region` | Ensembl gene ID |
| `phenotype` | Human-readable phenotype name (with `_int` suffix) |
| `phenocode` | Genebass/UKB phenocode |
| `mean_pheno_value` | SAIGE single-variant BETA, used as the APPV proxy |
| `SE` | Standard error of the BETA |
| `Pvalue` | SAIGE p-value |
| `AF` | Allele frequency in the cohort |
| `n_cases` | Number of cases/individuals with the phenotype |

### Notebook walkthrough

1. **Hail init** — Start Hail with GRCh38 as default reference.
2. **Load MT** — Read the Genebass MatrixTable and inspect schema. The MT is ~1 TB; no operation ever loads it fully into memory.
3. **Load associations** — Read the UKBBGym association parquet; extract unique Ensembl IDs and phenotype names.
4. **Gene ID mapping** — Map Ensembl IDs → gene symbols via the Biomart parquet; all 699 IDs map successfully.
5. **Column metadata** — Read phenotype metadata (description, phenocode, n_cases, …) directly from `MT_PATH/cols` as a standalone Hail Table (4,529 rows, fast) rather than calling `mt.cols()` which scans the full MT.
6. **Phenotype name normalisation** — Normalise both the MT `description` field and the UKBBGym phenotype names by stripping all non-alphanumeric characters and lowercasing, then match. 119/121 phenotypes are matched automatically; two impedance-based phenotypes (`weight_impedance_int`, `body_mass_index_bmi_impedance_int`) are unmatched because they lack a matching description in Genebass.
7. **Gene coordinates + symbols (GENCODE v29)** — Parse the v29 GTF for the 699 target genes' `(contig, start, end, gene_name)` into Polars `gene_coords`. v29 is the version Genebass annotated with (VEP v95), so its symbols match the MT's `gene` field and its boundaries match the variant assignments. (697/699 genes are in v29; 2 are newer IDs absent from v29 and from Genebass.)
8. **Extract entries (per phenotype, exact pairs)** — Loop over the 119 phenotypes. For each, take its associated genes, build padded `locus_interval`s from `gene_coords`, and: `filter_intervals` (locus-index prune) → `filter_rows` by **v29** gene symbol → `semi_join_cols` to that one phenotype → `.collect()` → Polars → one parquet shard per phenotype in `shards_allvars/`. Because each phenotype is paired only with its own genes, every collected row is a wanted `(gene, phenotype)` pair. ~119 small Hail jobs, no full scan, no pandas.
9. **Assemble & write parquet** — Fully lazy Polars: `scan_parquet(shards) → join n_cases → rename BETA→mean_pheno_value → sink_parquet` to the final output.

---

## Step 2 — Master table

**Notebook:** [utils/02_create_master_table.ipynb](utils/02_create_master_table.ipynb)

One table every analysis can use without a further join. Gene- and phenotype-level values (`phenotype`, `loftee_corr_dir`, …) are repeated on every variant row of the gene; that redundancy is deliberate, since it removes a join from every notebook. Only presentation metadata (`label`, `color`, `category`, `direction`) stays in [`../configs/`](../configs/).

**Grain:** one row per `(id, region)`.

**Scope:** genes FDR-significant in the regenie association file *and* carrying a non-NaN LOFTEE correlation, collapsed to the single strongest-|corr| phenotype per gene. Within those genes **no variant is filtered** — not by consequence, SNP/indel status, length, ClinVar, MAF or MAC. All of that is a downstream `.filter()`.

### Input data

| Variable | Path |
|---|---|
| `ANNO_PATH` | `<data_root>/variant_files/`<br>`qced_maf1e-3_loftee_olink_genes_EURunrelated/annotations_no_dup_20260630.parquet` — 248 annotation columns, produced by the main `../ukbb/` pipeline |
| `APPV_PATH` | `<data_root>/genebass/genebass_betas/genebass_betas_127phenos_allvars.parquet` — step 1 |
| `ASSOC_PATH` | `<data_root>/association_files/`<br>`regenie_127phenotypes_lofteeHC_mac20_EUR_miss20per.parquet` |
| `CORR_PATH` | `<data_root>/association_files/`<br>`regenie_127phenotypes_mac20_lofteeHC_EUR_correlations.parquet` |
| `CFG_DIR` | [`../configs/`](../configs/) |

Note the spelling of the data root on the Gagneur cluster: the mount holding the annotation file is spelled `ukbbgym` (two b's); the similarly named `ukbgym` mount is a different filesystem and does not hold it.

### Build steps

1. **`gene_trait_df`** — FDR ≤ 0.05 → inner join on the LOFTEE correlations → `drop_nans()` → keep the strongest-|corr| phenotype per gene. `loftee_corr_dir = corr / |corr|` carries the sign used for direction correction downstream; `loftee_corr_abs` picks the phenotype and is then dropped.
2. **`KEEP_COLS`** — derived from the YAMLs at runtime rather than hardcoded, so adding a tool to a config pulls its column in instead of silently producing nulls. Keep = every annotation named by any config (including the column names parsed out of `config_variant_classes.yaml` filter expressions) + the base columns that derived annotations need (`DERIVED`) + keys + all `consequence_*` flags + `amino_acids`/`protein_position` + the `*_is_na` masks (needed for the ClinVar notebook's imputation anti-join) + a short `EXTRA` list the notebooks use but no config names. Dropped: the four native `clinvar_*` flags — the category analyses recompute them from raw `clinical_significance`, and shipping both would mean two competing definitions — plus ~90 genuinely unused columns (`mirsvr-*`, `motif*`, `remapoverlap*`, `roulette-*`, `gerpn/gerps`, `grantham`, `aa_pos`, …). Result: **152 of 248** columns.
3. **Annotation slice** — `scan_parquet(ANNO_PATH).select(KEEP_COLS)` semi-joined to the retained genes.
4. **Join and write** — left join `gene_trait_df` on `region` (one phenotype per gene), then the betas on **`(id, phenotype)`**, adding `AC = AF × 2×394,841`, `AC_proxy = AF × 2 × n_cases`, `mean_pheno_value_dircor = mean_pheno_value × loftee_corr_dir` and the `is_cross_gene_beta` flag. Written with `sink_parquet(compression='zstd')`.

---

## Analysis notebooks

Every notebook under [`analysis/`](analysis/) reads the master table (`MASTER_PATH`, set at the
top of its first cell) and the shared [`../configs/`](../configs/) YAML — no
further joins against the annotation/association files. Some also read one or more of the
additional `data/` files in the inventory above; the *Needs* column says whether the published
master table alone is enough.

| Notebook | Produces | Needs |
|---|---|---|
| [`correlations.ipynb`](analysis/correlations.ipynb) | Pairwise missense-predictor heatmap: per-gene Spearman correlation with Genebass effect sizes, Wilcoxon signed-rank significance across gene–trait pairs | master table only |
| [`mean_phenotype.ipynb`](analysis/mean_phenotype.ipynb) | Mean direction-corrected carrier phenotype for missense variants stratified by protein structural / interaction annotation | master table only |
| [`protein_domains_correlations.ipynb`](analysis/protein_domains_correlations.ipynb) | Mean Spearman correlation within vs. outside structured domains, and the per-method paired difference, for TED domains and for AlphaFold2 pLDDT > 70 | master table only |
| [`noise_ceiling/noise_ceiling_one_region.ipynb`](analysis/noise_ceiling/noise_ceiling_one_region.ipynb) | Detectable-variance ceiling and variance captured per predictor, single unsplit population; noise estimated from AC=1 synonymous variants on the null-trait panel | + null-trait panel |
| [`noise_ceiling/noise_ceiling_TED_contrast.ipynb`](analysis/noise_ceiling/noise_ceiling_TED_contrast.ipynb) | Same ceiling framework, split within vs. outside TED domains | + null-trait panel |
| [`other_benchmarks/proteingym_snr.ipynb`](analysis/other_benchmarks/proteingym_snr.ipynb) | Signal-to-noise of pairwise predictor deltas, ProteinGym vs UKBBGym, percentile bootstrap | + ProteinGym |
| [`other_benchmarks/proteingym_correlations.ipynb`](analysis/other_benchmarks/proteingym_correlations.ipynb) | Per-method mean UKBBGym correlation (Genebass effect sizes) against mean ProteinGym correlation across human assays | + ProteinGym |
| [`other_benchmarks/clinvar_spearman_scatterplot.ipynb`](analysis/other_benchmarks/clinvar_spearman_scatterplot.ipynb) | Per-gene ClinVar pathogenicity auROC against per-gene Spearman correlation with the phenotype, one point per predictor | + ClinVar labels |
| [`other_benchmarks/expAssays_all_genes_correlations_pheno.ipynb`](analysis/other_benchmarks/expAssays_all_genes_correlations_pheno.ipynb) | Correlation with experimental deep mutational scanning assays (SGE, MaveDB) across genes | + all-variants master table, + ProteinGym, + DMS assays |

`proteingym_correlations.ipynb`'s UKBBGym-side correlation now reads directly off
the master table via the shared `gene_trait_tool_correlations` helper, like every other notebook
here, instead of separately re-selecting the best trait per gene from the raw association files
— the master table already applies that same selection at build time (see *Scope* in Step 2
above), so the gene-trait universe should be identical, but the published figure was produced by
the old code path and hasn't been re-run against this version to confirm bit-identical output.
