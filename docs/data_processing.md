# Data Processing

Before running the core analysis and benchmarking applets, the raw UK Biobank data must be processed into formatted Parquet files, and variant-level phenotype aggregates must be calculated.

## 1. Raw Data to Initial Parquets

The initial data processing steps—which convert UK Biobank raw BGEN/VCF formats and phenotype data into optimized Parquet formats—are maintained in the `ukb_gagneur` repository. 

> **Note:** The specifics of this upstream pipeline are abstracted here. The resulting output from this step is a set of initial Parquet files containing associations and basic annotations.

## 2. Average Phenotype Per Variant (APPV)

Once the upstream pipeline provides the initial association Parquet files, they must be converted into the specific format required by the UKBBGym analysis applets. This is done by computing the average phenotype value for each variant.

UKBBGym provides two dedicated DNAnexus applets for this step:
- `avg_pheno_per_variant_traits`: For quantitative trait phenotypes.
- `avg_pheno_per_variant_olink`: For Olink proteomics data.

### Running via the UKB-RAP GUI

1. Log into your UKB-RAP project.
2. Navigate to the **Applets** section or search for the `avg_pheno_per_variant` applets.
3. Select the applet and click **Run**.
4. In the configuration dialog, provide the required inputs (e.g., the initial Parquet files generated from the upstream pipeline).
5. Configure your desired parameters (such as minimum individuals per variant) and launch the job.

### Running via the Command Line Interface (CLI)

Alternatively, you can run the applet programmatically using the `dx-toolkit`.

**Example:**
```bash
dx run avg_pheno_per_variant_traits \
  -i input_parquet=project-XXX:/path/to/associations.parquet \
  -i min_individuals=20 \
  --destination project-XXX:/outputs/appv_traits/ \
  --brief
```

### Outputs

These applets generate the final `appv_parquet` files which contain:
- `id`: The variant identifier.
- `phenotype`: The specific trait or protein.
- `mean_pheno_value`: The aggregated average value.
- `n_individuals`: Number of individuals contributing to the aggregate.

These files serve as a core input for the [Analysis Applets](analysis.md).
