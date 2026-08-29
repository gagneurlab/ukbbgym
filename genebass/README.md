# Genebass Summary Statistics for UKBBGym

## Goal

Reconstruct the UKBBGym coding-variant benchmark from [Genebass](https://genebass.org/) single-variant summary statistics (394,841 UK Biobank exomes), so it can be run without individual-level data on the UK Biobank RAP. This is the setup that **runs locally**.

UKBBGym scores a variant by the **average phenotype among its carriers** (APPV) and evaluates variant effect predictors by the rank correlation between their scores and that quantity, per gene–trait association. Genebass publishes a SAIGE single-variant **BETA** instead. For a rare variant the two estimate the same thing up to a deterministic $(1-p)$ factor that cancels under any rank-based metric, so the BETA is used as the APPV proxy throughout.

Gene–trait pairs are *selected* by pLoF burden tests, but the benchmark scores variants of **every** consequence class, so all variants in the target genes are extracted, not just pLoF.

## Run this locally

You don't need the ~1 TB Genebass MatrixTable to reproduce the figures — the built master table
(`ukbbgym_genebass_20260824.parquet`) will be published on Hugging Face. Steps 1 and 2 below are
only needed if you want to rebuild it from scratch.

1. From the repo root: `uv sync` (see the root [README](../README.md) for environment details).
2. Download `ukbbgym_genebass_20260824.parquet` and place it at `data/ukbbgym_genebass_20260824.parquet`
   relative to the repo root — link to be added once the Hugging Face dataset is published.
   **This one file is enough for [`correlations_master_table.ipynb`](analysis/correlations_master_table.ipynb),
   [`mean_phenotype_master_table.ipynb`](analysis/mean_phenotype_master_table.ipynb) and
   [`protein_domains_correlations_master_file.ipynb`](analysis/protein_domains_correlations_master_file.ipynb).**
   Every other analysis notebook needs an additional file or two — see the table below and the
   *needs* column in *Analysis notebooks*.
3. Open any notebook under [`analysis/`](analysis/) with the `.venv` kernel and run top to bottom,
   or run all of them at once: `genebass/run_all.sh`. It executes every notebook under
   `analysis/` in place (via `nbconvert`, against the repo's `.venv`) and reports which ones
   failed — expected for any notebook past the first three until its extra `data/` inputs exist.
4. Figures are written to `paper_figures/` at the repo root.

### `data/` inventory

| File | Needed by | Hosting |
|---|---|---|
| `ukbbgym_genebass_20260824.parquet` | most `analysis/` notebooks (the master table) | Hugging Face (link TBD) |
| `ukbbgym_genebass_null_20260824.parquet` | `noise_ceiling/*` (unassociated/null-trait panel, AC=1 synonymous noise estimate) | not yet published |
| `ukbbgym_genebass_all_20260824.parquet` | `other_benchmarks/expAssays_all_genes_correlations_pheno.ipynb` (all-variants variant of the master table, not restricted to the 670 benchmark genes) | not yet published |
| `clinvar_significance_vep_annotations_processed_cadd_fill_na_20260804.parquet` | `other_benchmarks/clinvar_spearman_scatterplot_master_file.ipynb` (17,683-gene ClinVar label set; broader than the master table's `clinical_significance` column, see the notebook's own note) | not yet published |
| `proteingym_SNVs_with_readout_annotated_20260716.parquet` | `other_benchmarks/proteingym_snr_master_file.ipynb`, `other_benchmarks/proteingym_correlations_master_file.ipynb`, `other_benchmarks/expAssays_all_genes_correlations_pheno.ipynb` | not yet published |
| `DMS_Marsh_VEP.parquet`, `LDLR_Roth_Science_2025.parquet` | `other_benchmarks/expAssays_all_genes_correlations_pheno.ipynb` (additional experimental DMS assays) | not yet published |

Only the first row currently has a publication plan. The rest are flagged here rather than
silently left for a `FileNotFoundError` — treat any notebook past the first three as
not-yet-runnable until its inputs are published.

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
| 2 | [utils/02_create_master_table.ipynb](utils/02_create_master_table.ipynb) | `master_table_<date>.parquet` — one row per (variant, gene), everything joined; this is the file published as `ukbbgym_genebass_20260824.parquet` |

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
| [`correlations_master_table.ipynb`](analysis/correlations_master_table.ipynb) | Pairwise missense-predictor heatmap: per-gene Spearman correlation with Genebass effect sizes, Wilcoxon signed-rank significance across gene–trait pairs | master table only |
| [`mean_phenotype_master_table.ipynb`](analysis/mean_phenotype_master_table.ipynb) | Mean direction-corrected carrier phenotype for missense variants stratified by protein structural / interaction annotation | master table only |
| [`protein_domains_correlations_master_file.ipynb`](analysis/protein_domains_correlations_master_file.ipynb) | Mean Spearman correlation within vs. outside structured domains, and the per-method paired difference, for TED domains and for AlphaFold2 pLDDT > 70 | master table only |
| [`noise_ceiling/noise_ceiling_one_region.ipynb`](analysis/noise_ceiling/noise_ceiling_one_region.ipynb) | Detectable-variance ceiling and variance captured per predictor, single unsplit population; noise estimated from AC=1 synonymous variants on the null-trait panel | + null-trait panel |
| [`noise_ceiling/noise_ceiling_TED_contrast.ipynb`](analysis/noise_ceiling/noise_ceiling_TED_contrast.ipynb) | Same ceiling framework, split within vs. outside TED domains | + null-trait panel |
| [`other_benchmarks/proteingym_snr_master_file.ipynb`](analysis/other_benchmarks/proteingym_snr_master_file.ipynb) | Signal-to-noise of pairwise predictor deltas, ProteinGym vs UKBBGym, percentile bootstrap | + ProteinGym |
| [`other_benchmarks/proteingym_correlations_master_file.ipynb`](analysis/other_benchmarks/proteingym_correlations_master_file.ipynb) | Per-method mean UKBBGym correlation (Genebass effect sizes) against mean ProteinGym correlation across human assays | + ProteinGym |
| [`other_benchmarks/clinvar_spearman_scatterplot_master_file.ipynb`](analysis/other_benchmarks/clinvar_spearman_scatterplot_master_file.ipynb) | Per-gene ClinVar pathogenicity auROC against per-gene Spearman correlation with the phenotype, one point per predictor | + ClinVar labels |
| [`other_benchmarks/expAssays_all_genes_correlations_pheno.ipynb`](analysis/other_benchmarks/expAssays_all_genes_correlations_pheno.ipynb) | Correlation with experimental deep mutational scanning assays (SGE, MaveDB) across genes | + all-variants master table, + ProteinGym, + DMS assays |

`proteingym_correlations_master_file.ipynb`'s UKBBGym-side correlation now reads directly off
the master table via the shared `gene_trait_tool_correlations` helper, like every other notebook
here, instead of separately re-selecting the best trait per gene from the raw association files
— the master table already applies that same selection at build time (see *Scope* in Step 2
above), so the gene-trait universe should be identical, but the published figure was produced by
the old code path and hasn't been re-run against this version to confirm bit-identical output.
