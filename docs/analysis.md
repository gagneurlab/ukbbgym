# Benchmarking Applets

Once your data is processed, you can benchmark various variant annotation scoring models against LOFTEE correlations. This is handled by four dedicated DNAnexus applets.

## Overview of Applets

There are two primary types of analysis performed, each split into a version for traits and a version for Olink proteomics:

1. **Correlations Applets:** Compute how well each annotation score correlates with the average phenotype values for variants in a gene-trait association set.
    - `correlations_traits` (Quantitative traits)
    - `correlations_olink` (Olink proteomics)
2. **Top-N Variant Applets:** Rank variants by their annotation scores within each gene, compute cumulative phenotype effects at each rank, and compare across annotations.
    - `top_n_vars_traits`
    - `top_n_vars_olink`

## Required Inputs

All four applets generally require the following inputs:
- `associations_parquet`: Gene-trait/protein associations with LOFTEE correlation directions.
- `annotations_parquet`: Variant annotations including scores from all standard tools.
- `appv_parquet`: Average phenotype/protein per variant generated in the [Data Processing](data_processing.md) step.

## Running the Applets

### UKB-RAP GUI
1. Navigate to the Applet in your UKB-RAP project.
2. Select **Run** to open the configuration menu.
3. Attach the required Parquet files to the corresponding input fields.
4. Specify desired options (e.g., `variant_class`, `annotation_categories`, `mac_threshold`).
5. Run the job.

### CLI Example: `correlations_traits`

```bash
dx run correlations_traits \
  -i associations_parquet=project-XXX:/path/associations.parquet \
  -i annotations_parquet=project-XXX:/path/annotations.parquet \
  -i appv_parquet=project-XXX:/path/appv_traits.parquet \
  -i variant_class="missense" \
  -i annotation_categories="missense,conservation" \
  --destination project-XXX:/outputs/correlations/ \
  --brief
```

## Adding Custom Scores

A core feature of the analysis applets is the ability to benchmark your own custom scoring models against established tools.

You can upload custom model scores as a single Parquet file.
- **Required format:** Must contain an `id` column matching the variant IDs in your annotations and APPV files.
- **Score columns:** Any other float column is treated as a model score.

**CLI Example with Custom Scores:**
```bash
dx run top_n_vars_traits \
  -i associations_parquet=project-XXX:/path/associations.parquet \
  -i annotations_parquet=project-XXX:/path/annotations.parquet \
  -i appv_parquet=project-XXX:/path/appv_traits.parquet \
  -i custom_scores_parquet=project-XXX:/path/custom_scores.parquet \
  -i model_labels="MyCustomModel v1,MyCustomModel v2" \
  -i model_directions="1,-1" \
  --destination project-XXX:/outputs/top_n/
```

## Expected Outputs

The applets generate detailed metrics and publish-ready SVG plots:
- **results_parquet / cum_stats_parquet:** Tables containing all computed metrics, z-scores, and standard errors.
- **SVG Plots:** Boxplots, line plots, and pairwise Wilcoxon heatmaps for interpretation of the benchmarking results.
