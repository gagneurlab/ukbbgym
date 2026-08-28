# avg_pheno_per_variant_olink

Computes average phenotype per variant (appv) for Olink proteomics (PROTRIDER-corrected, no INT). Adds mean percentile rank by phenotype.

## Build

```bash
dx build applets/avg_pheno_per_variant_olink/ --destination <project>:/applets/ -f
```

## Inputs

| Input | Required | Type | Description |
|-------|----------|------|-------------|
| `associations_parquet` | Yes | file | Gene-protein associations (filtered by your criteria). Columns: `gene`. |
| `phenotypes_parquet` | Yes | file | PROTRIDER-corrected protein data. Columns: `sample`, protein gene IDs. |
| `annotations_parquet` | Yes | file | Variant annotations. Columns: `id`, `region`, `mac_ukb`. |
| `genotypes_parquet` | Yes | file | Genotypes in long format: `id`, `sample`, `gt`. |
| `samples_txt` | No | file | (Optional) Sample list to subset to (one sample ID per line, no header). If omitted, uses all samples from phenotypes file. |
| `mac_threshold` | No | int | (Optional) MAC filtering threshold. If omitted, no MAC filtering applied. |
| `pheno_chunk_size` | No | int | Chunk size for processing phenotypes. Default: `100`. |

## Output

`appv_parquet`: Parquet with columns `id`, `phenotype`, `n_individuals`, `mean_pheno_value`, `std_pheno_value`, `mean_pheno_value_rank`.

Filename: `average_pheno_per_variant_olink.parquet`

- `mean_pheno_value_rank`: descending rank (1 = highest mean value)

## Details

- Phenotypes are PROTRIDER-corrected and already normalized (no INT applied).
- All associations from the provided file are used (no p-value filtering).
- Instance: `mem3_ssd1_v2_x16` (128 GB). Runtime ~20 minutes.
