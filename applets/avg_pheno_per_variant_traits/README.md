# avg_pheno_per_variant

Computes average phenotype per variant (appv) for quantitative traits with inverse-normal transformation (INT using Blom's constant).

## Build

```bash
dx build applets/avg_pheno_per_variant/ --destination <project>:/applets/ -f
```

## Inputs

| Input | Required | Type | Description |
|-------|----------|------|-------------|
| `associations_parquet` | Yes | file | Gene-trait associations (filtered by your criteria). Columns: `region`, `phenotype`. |
| `phenotypes_parquet` | Yes | file | Phenotype data. Columns: `individual`, all phenotype columns. |
| `annotations_parquet` | Yes | file | Variant annotations. Columns: `id`, `region`, `mac_ukb`. |
| `genotypes_parquet` | Yes | file | Genotypes in long format: `id`, `sample`, `gt`. |
| `samples_txt` | No | file | (Optional) Sample list to subset to (one sample ID per line, no header). If omitted, uses all samples from phenotypes file. |
| `mac_threshold` | No | int | (Optional) MAC filtering threshold. If omitted, no MAC filtering applied. |
| `pheno_chunk_size` | No | int | Chunk size for processing phenotypes. Default: `10`. |

## Output

`appv_parquet`: Parquet with columns `id`, `phenotype`, `n_individuals`, `mean_pheno_value`, `std_pheno_value`.

Filename: `average_pheno_per_variant_traits.parquet`

## Details

- INT uses Blom's constant (3/8) for rank normalization.
- All associations from the provided file are used (no FDR filtering).
- Instance: `mem3_ssd1_v2_x16` (128 GB). Runtime ~5 minutes.
