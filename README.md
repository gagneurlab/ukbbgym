# UKBBGym

<img src="assets/logo.png" alt="UKBBGym" width="140" align="right">

Population-scale phenotype benchmark for variant effect predictors.

UKBBGym evaulates variant scoring methods against **phenotypes observed in human
carriers**, rather than against curated clinical labels or cell-culture assays. For every
ultra-rare variant (allele count ≤ 20) within ±5 kb of a gene body, the phenotype is averaged
across its heterozygous carriers, after correction for covariates and common-variant polygenic
risk, and predictors are scored by the rank correlation between their scores and these
per-variant carrier means. Two phenotype tracks: **670 gene–trait associations** across 95
quantitative traits, and **1,076 gene–protein abundance associations** from Olink plasma
proteomics.

## Which setup do I want?

The benchmark can be run three ways, depending on what data you have access to:

| Setup | Data source | Runs locally | Access |
|---|---|---|---|
| [`genebass/`](genebass/README.md) | [Genebass](https://genebass.org/) single-variant summary stats (394,841 UKB exomes) | **Yes** | master table published on Hugging Face |
| [`all_x_all/`](all_x_all/README.md) | All of Us All-by-All summary stats | No | [All of Us Researcher Workbench](https://www.researchallofus.org/) only |
| [`ukbb/`](ukbb/README.md) | Individual-level UK Biobank (the primary benchmark) | No | [UK Biobank RAP](https://ukbiobank.dnanexus.com/) only |

`genebass/` is the only setup you can clone and run on a laptop — start there unless you
specifically have UKBB RAP or All of Us Workbench access.
To run the benchmark on the UKBB RAP, genotype, Olink, and phenotype data processing can be done using applets from this git repository: [gagneurlab/ukbb-utils-gagneur](https://github.com/gagneurlab/ukbb-utils-gagneur)

## Start here (genebass, local)

```bash
git clone <this-repo-url>
cd ukbbgym
uv sync                       # creates .venv/ from uv.lock — see "Environment setup" below

# Log in to Hugging Face (the master table lives in the private gagneurlab/ukbbgym dataset):
.venv/bin/huggingface-cli login          # or: export HF_TOKEN=...

# Run any of these three with the .venv kernel — no manual download needed, they fetch
# the master table into data/ (gitignored) on first run:
#   genebass/analysis/correlations.ipynb
#   genebass/analysis/mean_phenotype.ipynb
#   genebass/analysis/protein_domains_correlations.ipynb
# Every other genebass/analysis notebook needs one or more additional input files —
# see the data/ inventory in genebass/README.md.
```

See [`genebass/README.md`](genebass/README.md) for the method, the notebook-by-figure table,
and how to rebuild the master table from raw Genebass/UKB inputs instead of the published one.

## Repository layout

```
genebass/       Genebass summary-statistics reconstruction — runs locally (see above)
all_x_all/      All of Us All-by-All reconstruction — correlation heatmap only, Workbench-only
ukbb/           primary UK Biobank RAP pipeline — applets, docs, RAP-only
configs/        YAML shared by genebass/+all_x_all, plus configs/ukbb/ for the RAP pipeline —
                predictor sets, variant-class filters, colors/labels
utils/          code shared across setups: variant_filtering.py, REGENIE burden testing,
                annotation regeneration, experimental-assay data
envs/           Hail/Spark environment spec (separate from the core analysis environment)
```

Each directory has its own README with the detail specific to it.

## Environment setup

### Core analysis environment (all figure notebooks)

The default install **assumes the master table already exists** — it covers reading that table
and reproducing the figures, not producing it. With [uv](https://docs.astral.sh/uv/)
(`curl -LsSf https://astral.sh/uv/install.sh | sh`):

```bash
uv sync                       # creates .venv/ from uv.lock, exact pinned versions
```

When running the notebooks, pick `.venv/bin/python` as the kernel.

**Without the lockfile.** [`pyproject.toml`](pyproject.toml) is self-contained; `uv.lock` only
records the exact versions that resolution picked, so the figures can be reproduced bit-for-bit
later. To resolve fresh against the declared ranges instead — with uv, or with plain pip and no
uv at all:

```bash
uv pip install -e .                # into the active environment; add ".[annotation]" for the extra
```

**Conda/mamba equivalent:**
```bash
mamba env create -f envs/environment.yml
```

### Optional: `annotation` extra (regenerating variant annotations)

| Extra | Install | Needed by |
|---|---|---|
| `annotation` | `uv sync --extra annotation` | [`utils/annotations/*`](utils/README.md) — adds `pyBigWig`, `pyfaidx`, `pysam`, `biopython` for reading bigWig / FASTA / tabix reference files, plus `duckdb` and `requests` for gene subsetting and record fetching |

### Optional: Hail environment (summary-statistics extraction)

Needed only to re-extract the summary statistics that feed the master table — not needed if
you're using the published Hugging Face table. The Hail/Spark notebooks —
[`genebass/utils/01_read_hail_sumstats.ipynb`](genebass/utils/01_read_hail_sumstats.ipynb) and
[`all_x_all/utils/*`](all_x_all/README.md) — need a **separate** environment. Hail pins
`pyspark` and an older `numpy` range that would drag the analysis stack backwards, so it's kept
out of `pyproject.toml`:

```bash
uv venv --python 3.11 .venv-hail
uv pip install --python .venv-hail -r envs/requirements-hail.txt
```

Requires a **Java 11** JDK on `PATH` (`java -version` to check; Temurin 11 is what Hail
recommends) — Hail 0.2.139 runs on Spark 3.5 and will not start against a newer JDK. These
notebooks read TB-scale MatrixTables and are normally run on a cluster with the data mounted
(or the All of Us Workbench, for `all_x_all/`), not on a laptop.

## Data access

- **genebass** — the built master table (`genebass_annotated.parquet`) is published in the
  private [`gagneurlab/ukbbgym`](https://huggingface.co/datasets/gagneurlab/ukbbgym) Hugging
  Face dataset and fetched automatically by the notebooks (login required, see below); the
  Genebass MatrixTable it can optionally be rebuilt from is available from the Genebass authors
  on Google Cloud. See [`genebass/README.md`](genebass/README.md).
- **all_x_all** — inputs live only inside the All of Us Researcher Workbench and cannot be
  downloaded; see [`all_x_all/README.md`](all_x_all/README.md).
- **ukbb** — inputs are individual-level UK Biobank data, accessible only on the RAP under your
  own UKB application; see [`ukbb/README.md`](ukbb/README.md).

## Contributing

Install pre-commit before committing changes, so notebook output stripping and formatting run
automatically:

```bash
pip install pre-commit      # or: sudo apt install pre-commit
pre-commit install
```
