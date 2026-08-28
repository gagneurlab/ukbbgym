#!/usr/bin/env python3
"""
add_missense_variant_annotations.py
-----------------------------------
Annotate exome/missense variants with effect-prediction scores.

Ported from the DNAnexus applet to run standalone on the All of Us Researcher
Workbench (no dxpy; input/output are plain filesystem or gs:// paths).

Takes the output of vep_loftee_parallel (VEP + LOFTEE + gnomAD AFs + structural features) and adds:

  Step 1 — Downloaded at runtime:
    am_pathogenicity   (AlphaMissense, Zenodo ~1.2 GB)
    clinvar / clinvar_patho dummies  (ClinVar NCBI FTP ~600 MB)

  Step 3 — Auto-downloaded scores (all attempted; skipped gracefully on failure):
    revel_score (REVEL v1.3 website)
    clinpred_score (ClinPred hg38)
    bayes_del (BayesDel noAF)
    cpt1_llr (CPT-1 per-protein, Zenodo 8140323)
    gpn_score (GPN-MSA, remote tabix via HuggingFace bgz)
    cadd_raw (CADD v1.7 GRCh38, remote tabix; SNVs from whole_genome_SNVs.tsv.gz,
              indels from the gnomAD v4.0 genomes indel file)
    phylop_100way (UCSC bigWig, remote pyBigWig)

  Step 4 — Derived columns:
    indel length categories (1bp_del, 1bp_ins, 2_5bp_del, ...),

  Step 5 — Fill nulls + _is_na columns

Note: VEP-derived structural features (relative_cds_position, indel flags, dist_to_tss,
gene_length, gene_name, next_in_frame_relative) are now computed in vep_loftee_parallel applet.
Protein domain annotations (MobiDB, TED, low-complexity) were removed to avoid user file inputs.

Logic ported from:
  deeprvat_wgs/annotation/annotation_functions.py
  ukbgym/utils/annotations/add_ukbgym_annotations_to_bcf2parquet.ipynb
"""

import gc
import gzip
import logging
import os
import shutil
import subprocess
import sys
import tarfile
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import polars as pl
import yaml

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def _find_tabix() -> str:
    found = shutil.which("tabix")
    if found:
        return found
    raise FileNotFoundError(
        "tabix not found on PATH. "
        "Install via: pip install pysam  (provides tabix) "
        "or conda install -c bioconda htslib"
    )


_TABIX: str | None = None


def _tabix_bin() -> str:
    global _TABIX
    if _TABIX is None:
        _TABIX = _find_tabix()
        logger.info(f"Using tabix: {_TABIX}")
    return _TABIX


def _resolve_work_dir() -> str:
    """
    Pick a writable directory for caching downloads.
    Priority: $ANNO_WORK_DIR env var → ~/anno_work (AoU Jupyter home is
    persistent and roomy) → /tmp/anno_work.
    """
    override = os.environ.get("ANNO_WORK_DIR")
    if override:
        return override
    candidates = [
        os.path.expanduser("~/anno_work"),
        "/tmp/anno_work",
    ]
    for c in candidates:
        try:
            os.makedirs(c, exist_ok=True)
            test = os.path.join(c, ".write_test")
            with open(test, "w") as fh:
                fh.write("")
            os.remove(test)
            return c
        except (OSError, PermissionError):
            continue
    return "/tmp"


WORK_DIR = _resolve_work_dir()

# ── Download URLs ──────────────────────────────────────────────────────────
AM_URL = "https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz"
CLINVAR_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
)

REVEL_URL = (
    # "https://rothsupport.hms.harvard.edu/revel-v1.3_all_chromosomes.zip"
    "https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip?download=1"
)
# ClinPred hg38 precomputed scores — distributed via Google Drive
# Source: https://drive.google.com/file/d/1e0kd9hO1uCEuGzAhwEmLCNluqDFJE28y/
CLINPRED_GDRIVE_ID = "1e0kd9hO1uCEuGzAhwEmLCNluqDFJE28y"
# BayesDel noAF whole-genome precomputed scores — distributed via Google Drive
# (the original fenglab.chpc.utah.edu URL has been taken offline).
# Source link: https://drive.google.com/file/d/1TcRxedFOuLwDvenWxo4pPVN8STqy-s7W/
BAYESDEL_GDRIVE_ID = "1TcRxedFOuLwDvenWxo4pPVN8STqy-s7W"
# UniProt human ID mapping file
UNIPROT_IDMAP_URL = (
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
    "knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz"
)
# popEVE / EVE / ESM1v per-variant VCF (GRCh38, ~1.4 GB)
POPEVE_URL = (
    "https://data.evemodel.org/popeve/v1.1/downloads/grch38_popEVE_ukbb_20250715.vcf.gz"
)
# CPT-1 per-protein scores from Zenodo (record 8140323)
CPT1_ZENODO_URLS = [
    "https://zenodo.org/records/8140323/files/CPT1_score_EVE_set.zip?download=1",
    "https://zenodo.org/records/8140323/files/CPT1_score_no_EVE_set_1.zip?download=1",
    "https://zenodo.org/records/8140323/files/CPT1_score_no_EVE_set_2.zip?download=1",
]
# GPN-MSA bgzip+tabix TSV from HuggingFace
GPN_MSA_URL = (
    "https://huggingface.co/datasets/songlab/gpn-msa-hg38-scores/resolve/main/"
    "scores.tsv.bgz"
)
GPN_MSA_TBI_URL = GPN_MSA_URL + ".tbi"

# CADD GRCh38 precomputed SNV and indel scores (remote tabix)
CADD_SNV_URL = (
    "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/"
    "whole_genome_SNVs.tsv.gz"
)
# CADD has no whole-genome indel score set; the closest coverage is this
# gnomAD v4.0 genomes indel file.
CADD_INDEL_URL = (
    "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/"
    "gnomad.genomes.r4.0.indel.tsv.gz"
)
# PhyloP 100-way vertebrate conservation (UCSC bigWig, accessed remotely)
PHYLOP_100WAY_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/"
    "hg38.phyloP100way.bw"
)
PHYLOP_17WAY_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP17way/"
    "hg38.phyloP17way.bw"
)
PHYLOP_30WAY_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP30way/"
    "hg38.phyloP30way.bw"
)
# AlphaFold human proteome v6 — per-protein PDB.gz tarball (~4.8 GB)
# Each archive entry is named AF-{uniprot_id}-F1-model_v6.pdb.gz
ALPHAFOLD_HUMAN_TAR_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/"
    "UP000005640_9606_HUMAN_v6.tar"
)
# PIONEER High-confidence interface residues (Xiong et al., Nat Biotech 2024)
PIONEER_HIGH_URL = "https://pioneer.yulab.org/static/predictions/high/human.txt"
# Gencode v40 human annotation GTF (~50 MB compressed) — used by REVEL to build
# a complete ENST → ENSG transcript→gene map (Gencode covers every transcript,
# unlike UniProt idmapping which only lists canonical transcripts per protein).
GENCODE_GTF_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/"
    "gencode.v40.annotation.gtf.gz"
)
GENCODE_FASTA_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/"
    "GRCh38.primary_assembly.genome.fa.gz"
)

# ── Utilities ──────────────────────────────────────────────────────────────
def run(cmd: str) -> None:
    logger.info(f"$ {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def aria2c_download(url: str, outpath: str, max_retries: int = 5) -> bool:
    """Download a file with aria2c, resuming on TLS resets or connection drops.

    The CADD server (krishna.gs.washington.edu) and others drop long-running
    TLS connections mid-transfer. This wrapper retries up to `max_retries`
    times, always passing --continue=true so aria2c resumes from where it left
    off rather than restarting.

    Servers that don't support HTTP range requests (e.g. pioneer.yulab.org)
    are detected via a known-host list and downloaded with a single connection
    to avoid noisy per-segment 416/range errors.

    Returns True only when the file exists AND has no .aria2 control file
    (which signals an incomplete transfer).
    """
    # Servers known to ignore or reject range requests — use single connection.
    no_range_hosts = ("pioneer.yulab.org",)
    if any(h in url for h in no_range_hosts):
        connections = "-x 1 -s 1"
    else:
        connections = "-x 16 -s 16"

    os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
    outdir = os.path.dirname(outpath) or "."
    outfile = os.path.basename(outpath)
    aria2_ctrl = outpath + ".aria2"

    for attempt in range(1, max_retries + 1):
        # Already complete from a previous run
        if os.path.exists(outpath) and not os.path.exists(aria2_ctrl):
            logger.info(f"  Already complete: {outpath}")
            return True

        if attempt > 1:
            size_mb = os.path.getsize(outpath) / 1024**2 if os.path.exists(outpath) else 0
            logger.info(f"  Retry {attempt}/{max_retries} — resuming from {size_mb:.0f} MB")
        try:
            run(
                f"aria2c {connections} "
                f"--continue=true "
                f"--allow-overwrite=true "
                f"--max-tries=3 "
                f"--retry-wait=5 "
                f"--timeout=60 "
                f"--connect-timeout=30 "
                f"'{url}' -d '{outdir}' -o '{outfile}'"
            )
        except subprocess.CalledProcessError:
            if attempt == max_retries:
                logger.warning(f"  Download failed after {max_retries} attempts: {url}")
                return False
            logger.info(f"  aria2c exited with error; will retry ({attempt}/{max_retries})")
            time.sleep(10 * attempt)
            continue

        if os.path.exists(outpath) and not os.path.exists(aria2_ctrl):
            return True
        logger.info(f"  Transfer incomplete after attempt {attempt}; retrying...")

    return False

def _tabix_download(url: str, outpath: str) -> bool:
    """Choose the right download strategy for a tabix bgz file.

    HuggingFace (and other CDNs that use pre-signed S3 URLs with per-chunk
    ByteRange policies) are incompatible with aria2c's multi-connection
    splitting: each parallel segment gets a signed byte range baked into its
    policy, so when aria2c tries to resume with a different offset it gets 403.

    For these URLs we fall back to a single-connection curl, which follows
    redirects, supports byte-range resuming with -C -, and doesn't need the
    content-length up front. For everything else we use aria2c as normal.
    """
    huggingface_hosts = ("huggingface.co", "hf.co", "xethub.hf.co")
    is_hf = any(h in url for h in huggingface_hosts)

    if is_hf:
        logger.info("  Using curl (single-connection) for HuggingFace pre-signed URL")
        aria2_ctrl = outpath + ".aria2"
        # Clean up any stale aria2 control file that would confuse the
        # completion check but doesn't help curl.
        if os.path.exists(aria2_ctrl):
            os.remove(aria2_ctrl)
        for attempt in range(1, 6):
            if os.path.exists(outpath):
                size_mb = os.path.getsize(outpath) / 1024**2
                logger.info(f"  curl attempt {attempt}: resuming from {size_mb:.0f} MB")
            try:
                run(
                    f"curl -L --retry 5 --retry-delay 10 --retry-connrefused "
                    f"-C - "           # resume from existing partial file
                    f"--max-time 7200 "  # 2 h ceiling; 37 GB at ~5 MB/s ≈ 2 h
                    f"-o '{outpath}' "
                    f"'{url}'"
                )
                # curl exits 0 even on partial success with -C -; verify
                # the file grew and no .aria2 control file exists
                if os.path.exists(outpath) and os.path.getsize(outpath) > 0:
                    return True
            except subprocess.CalledProcessError:
                if attempt == 5:
                    logger.warning(f"  curl failed after 5 attempts: {url}")
                    return False
                logger.info(f"  curl attempt {attempt} failed; retrying in 30s")
                time.sleep(30)
        return False
    else:
        return aria2c_download(url, outpath)


def _ensure_local_tabix(url: str, local_dir: str | None = None) -> str | None:
    """
    Ensure a tabix-indexed bgz file is available locally (data + .tbi).
    Downloads via aria2c if not already cached. Returns the local data path,
    or None if download failed.
    """
    if local_dir is None:
        local_dir = WORK_DIR
    basename = os.path.basename(url.split("?")[0])
    local_path = os.path.join(local_dir, basename)
    tbi_local = local_path + ".tbi"
    tbi_url = url + ".tbi"

    # aria2c writes a <file>.aria2 control file while downloading; it removes
    # it only on successful completion. If it exists the file is incomplete.
    aria2_ctrl = local_path + ".aria2"
    need_download = not os.path.exists(local_path) or os.path.exists(aria2_ctrl)

    if need_download:
        size_str = ""
        try:
            import urllib.request  # noqa: PLC0415
            with urllib.request.urlopen(url, timeout=15) as resp:
                cl = resp.headers.get("Content-Length")
                if cl:
                    size_str = f" ({int(cl) / 1e9:.1f} GB)"
        except Exception:
            pass
        if os.path.exists(aria2_ctrl):
            logger.info(
                f"  Resuming incomplete download{size_str}: {local_path}"
            )
        else:
            logger.info(f"  Downloading tabix data{size_str}: {url}")
        t0 = time.time()
        if not _tabix_download(url, local_path):
            logger.warning(f"  tabix data download failed: {url}")
            return None
        logger.info(f"  Downloaded {basename} in {time.time() - t0:.0f}s")
    else:
        logger.info(f"  Using cached tabix data: {local_path}")

    if not os.path.exists(tbi_local):
        logger.info(f"  Downloading tabix index: {tbi_url}")
        if not aria2c_download(tbi_url, tbi_local):
            logger.warning(f"  tabix index download failed: {tbi_url}")
            return None

    return local_path


def gdrive_download(file_id: str, outpath: str) -> bool:
    """Download a file from Google Drive given its file ID.

    Handles the large-file virus-scan confirmation by parsing the form Google
    serves on the first hit and re-issuing the request to
    `drive.usercontent.google.com/download`.  Returns True on success.
    """
    import re
    from urllib.parse import urlparse, parse_qs

    # Already present from a previous run — skip the network round-trip.
    # 1024 mirrors the "too small / error page" guard below.
    if os.path.exists(outpath) and os.path.getsize(outpath) > 1024:
        logger.info(f"  Already present: {outpath}")
        return True

    try:
        import requests
    except ImportError:
        logger.warning("Google Drive download needs the `requests` package")
        return False

    os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
    url = f"https://drive.google.com/uc?id={file_id}&export=download"
    session = requests.Session()
    CONNECT_TIMEOUT = 60  # seconds

    try:
        logger.info(f"Google Drive: fetching confirmation page for {file_id}")
        r = session.get(url, allow_redirects=True, timeout=CONNECT_TIMEOUT)
        ctype = r.headers.get("content-type", "").lower()
        logger.info(f"Google Drive: initial response content-type={ctype!r} status={r.status_code}")

        if "text/html" in ctype:
            text = r.text
            m = re.search(
                r'action="(https://drive\.usercontent\.google\.com/download[^"]*)"',
                text,
            )
            if not m:
                logger.warning(
                    f"Google Drive: could not find download form for {file_id}. "
                    f"Page snippet: {text[:500]!r}"
                )
                return False
            action = m.group(1).replace("&amp;", "&")
            params = {}
            for inp in re.finditer(
                r'<input[^>]+name="([^"]+)"[^>]+value="([^"]*)"', text
            ):
                params[inp.group(1)] = inp.group(2)
            parsed = urlparse(action)
            base = f"{parsed.scheme}://{parsed.netloc}{parsed.path}"
            merged = {k: v[0] for k, v in parse_qs(parsed.query).items()}
            merged.update(params)
            logger.info(f"Google Drive: submitting form to {base} with params {list(merged.keys())}")
            r = session.get(base, params=merged, stream=True, allow_redirects=True, timeout=CONNECT_TIMEOUT)
            logger.info(f"Google Drive: download response status={r.status_code} content-type={r.headers.get('content-type')!r} size={r.headers.get('content-length', 'unknown')} bytes")
        else:
            r = session.get(url, stream=True, allow_redirects=True, timeout=CONNECT_TIMEOUT)

        chunk_size = 1 << 20  # 1 MB
        written = 0
        with open(outpath, "wb") as f:
            for chunk in r.iter_content(chunk_size):
                if chunk:
                    f.write(chunk)
                    written += len(chunk)
                    if written % (50 << 20) < chunk_size:
                        logger.info(f"Google Drive: downloaded {written >> 20} MB so far ...")
        logger.info(f"Google Drive: finished writing {written >> 20} MB to {outpath}")
    except Exception as e:
        logger.warning(f"Google Drive download failed for {file_id}: {e}")
        return False

    size = os.path.getsize(outpath) if os.path.exists(outpath) else 0
    if size < 1024:
        logger.warning(f"Google Drive download too small ({size} bytes) for {file_id}; likely an error page")
        return False
    return True


def unzip(zip_path: str, extract_dir: str) -> str:
    """Extract a .zip, .tar.gz, or plain .gz archive; return the extraction directory.

    Idempotent: writes a per-archive sentinel ({extract_dir}/.extracted_<basename>)
    after a successful extraction and skips re-extraction if it already exists. The
    marker is written only on completion, so a partial extraction is never mistaken
    for a complete one.
    """
    os.makedirs(extract_dir, exist_ok=True)
    marker = os.path.join(extract_dir, f".extracted_{os.path.basename(zip_path)}")
    if os.path.exists(marker):
        logger.info(f"  Already extracted: {zip_path} -> {extract_dir}")
        return extract_dir

    logger.info(f"Extracting {zip_path} to {extract_dir}")
    with open(zip_path, "rb") as fh:
        magic = fh.read(4)
    if magic[:2] == b"\x1f\x8b":
        # gzip — could be .tar.gz or a plain compressed file
        if tarfile.is_tarfile(zip_path):
            with tarfile.open(zip_path, "r:gz") as tf:
                tf.extractall(extract_dir)
        else:
            # plain .gz — decompress single file
            out_name = os.path.splitext(os.path.basename(zip_path))[0]
            with gzip.open(zip_path, "rb") as gz_in, open(os.path.join(extract_dir, out_name), "wb") as f_out:
                shutil.copyfileobj(gz_in, f_out)
    else:
        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(extract_dir)
    # Mark complete so reruns skip re-extraction.
    with open(marker, "w") as fh:
        fh.write("")
    return extract_dir


def _idempotent_join(
    left: pl.LazyFrame,
    right: pl.LazyFrame,
    on: list[str],
    how: str = "left",
    **kwargs,
) -> pl.LazyFrame:
    """Left-join that is safe to re-run in a notebook.

    Drops any columns the right frame would add that already exist in the left
    frame (excluding the join keys), so repeated cell executions don't produce
    DuplicateError.
    """
    left_cols = set(left.collect_schema().names())
    right_cols = set(right.collect_schema().names())
    keys = set(on) if isinstance(on, list) else {on}
    to_drop = (right_cols - keys) & left_cols
    if to_drop:
        left = left.drop(list(to_drop))
    return left.join(right, on=on, how=how, **kwargs)


def _build_uniprot_map(work_dir: str) -> pl.DataFrame:
    """
    Download and parse the UniProt human idmapping file.
    Returns a DataFrame with columns:
      region     (Ensembl gene ID, e.g. ENSG00000141510)
      uniprot_id (UniProt accession, e.g. P04637)
      gene_name  (HGNC gene symbol, e.g. TP53)
      entry_name (UniProt Entry Name / mnemonic, e.g. P53_HUMAN)
    Used by CPT-1 and protein-domain steps.
    """
    path = os.path.join(work_dir, "HUMAN_9606_idmapping.dat.gz")
    ok = aria2c_download(UNIPROT_IDMAP_URL, path)
    if not ok:
        logger.warning("UniProt idmapping download failed; returning empty map")
        return pl.DataFrame(
            schema={"region": pl.Utf8, "uniprot_id": pl.Utf8, "gene_name": pl.Utf8, "entry_name": pl.Utf8}
        )

    df = pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        new_columns=["uniprot_id", "id_type", "value"],
        schema_overrides={"uniprot_id": pl.Utf8, "id_type": pl.Utf8, "value": pl.Utf8},
    )

    ensembl_map = (
        df.filter(pl.col("id_type") == "Ensembl")
        .select(["uniprot_id", "value"])
        .with_columns(region=pl.col("value").str.split(".").list.first())
        .select(["uniprot_id", "region"])
    )
    gene_name_map = (
        df.filter(pl.col("id_type") == "Gene_Name")
        .select(["uniprot_id", "value"])
        .rename({"value": "gene_name"})
        .group_by("uniprot_id")
        .first()
    )
    entry_name_map = (
        df.filter(pl.col("id_type") == "UniProtKB-ID")
        .select(["uniprot_id", "value"])
        .rename({"value": "entry_name"})
        .group_by("uniprot_id")
        .first()
    )
    return (
        ensembl_map
        .join(gene_name_map, on="uniprot_id", how="left")
        .join(entry_name_map, on="uniprot_id", how="left")
        .unique(subset=["region"])
    )


# ── Step 2: Downloaded annotations ────────────────────────────────────────


def _build_am_id_maps(idmap_path: str) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Build the two ID maps used by the AlphaMissense merge.

    Returns:
      uniprot_to_ensg : (uniprot_id, region) — FULL map, no dedup. Captures every
                        UniProt accession listed in idmapping with an Ensembl gene.
      enst_to_ensg    : (enst_id, region)    — derived via UniProt as bridge,
                        version suffix stripped.
    Both maps are needed because AM `uniprot_id`s and `transcript_id`s have
    partially overlapping coverage of the human proteome.
    """
    df = pl.read_csv(
        idmap_path,
        separator="\t",
        has_header=False,
        new_columns=["uniprot_id", "id_type", "value"],
        schema_overrides={"uniprot_id": pl.Utf8, "id_type": pl.Utf8, "value": pl.Utf8},
    )
    uniprot_to_ensg = (
        df.filter(pl.col("id_type") == "Ensembl")
        .select(["uniprot_id", "value"])
        .with_columns(region=pl.col("value").str.split(".").list.first())
        .select(["uniprot_id", "region"])
        .unique()
    )
    enst_uid = (
        df.filter(pl.col("id_type") == "Ensembl_TRS")
        .select(["uniprot_id", "value"])
        .with_columns(enst_id=pl.col("value").str.split(".").list.first())
        .select(["uniprot_id", "enst_id"])
    )
    enst_to_ensg = (
        enst_uid.join(uniprot_to_ensg, on="uniprot_id", how="inner")
        .select(["enst_id", "region"]).unique()
    )
    return uniprot_to_ensg, enst_to_ensg

# ── Group 1: Missense / protein-function predictors ──


def step1a_alphamissense(
    annos: pl.LazyFrame, am_tsv_path: str, work_dir: str
) -> pl.LazyFrame:
    """Add AlphaMissense pathogenicity scores.

    Joins on (chrom, pos, ref, alt, region) where region is the Ensembl gene ID.
    Resolves AM rows to ENSG via TWO independent maps and coalesces:
      1. Primary  : AM `uniprot_id`     → ENSG (full uniprot→ENSG, no dedup)
      2. Fallback : AM `transcript_id`  → ENSG (via Ensembl_TRS bridge)
    The full uniprot map captures more genes than the deduped version used
    elsewhere; the transcript map fills any gaps where AM's uniprot_id is
    missing from idmapping (e.g., deprecated accessions).
    """
    logger.info("Step 1a: Merging AlphaMissense scores")
    idmap_path = os.path.join(work_dir, "HUMAN_9606_idmapping.dat.gz")
    uniprot_to_ensg, enst_to_ensg = _build_am_id_maps(idmap_path)

    base = pl.scan_csv(
        am_tsv_path, separator="\t", skip_rows=3,
        schema_overrides={
            "#CHROM": pl.Utf8, "POS": pl.Int64, "am_pathogenicity": pl.Float32,
        },
    ).rename({"#CHROM": "chrom", "POS": "pos", "REF": "ref", "ALT": "alt"})

    am_via_uniprot = (
        base.select(["chrom", "pos", "ref", "alt", "uniprot_id", "am_pathogenicity"])
        .join(uniprot_to_ensg.lazy(), on="uniprot_id", how="left")
        .filter(pl.col("region").is_not_null())
        .group_by(["chrom", "pos", "ref", "alt", "region"])
        .agg(pl.col("am_pathogenicity").max().alias("am_uniprot"))
    )
    am_via_transcript = (
        base.select(["chrom", "pos", "ref", "alt", "transcript_id", "am_pathogenicity"])
        .with_columns(enst_id=pl.col("transcript_id").str.split(".").list.first())
        .join(enst_to_ensg.lazy(), on="enst_id", how="left")
        .filter(pl.col("region").is_not_null())
        .group_by(["chrom", "pos", "ref", "alt", "region"])
        .agg(pl.col("am_pathogenicity").max().alias("am_transcript"))
    )
    annos = (
        annos
        .join(am_via_uniprot,    on=["chrom", "pos", "ref", "alt", "region"], how="left")
        .join(am_via_transcript, on=["chrom", "pos", "ref", "alt", "region"], how="left")
        .with_columns(am_pathogenicity=pl.coalesce(["am_uniprot", "am_transcript"]))
        .drop(["am_uniprot", "am_transcript"])
    )
    logger.info("  AlphaMissense OK")
    return annos


def _build_enst_to_ensg_from_gtf(gtf_path: str) -> pl.DataFrame:
    """Build ENST → ENSG map from a Gencode GTF.

    Gencode covers every transcript (canonical + alternative), so this map is
    far more complete than the UniProt idmapping bridge used elsewhere.
    Strips version suffixes from both `transcript_id` and `gene_id`.
    Returns columns: `enst_id`, `region`.
    """
    return (
        pl.read_csv(
            gtf_path,
            separator="\t",
            comment_prefix="#",
            has_header=False,
            new_columns=[
                "chrom", "source", "feature", "start", "end", "score",
                "strand", "frame", "attributes",
            ],
            schema_overrides={"feature": pl.Utf8, "attributes": pl.Utf8},
            ignore_errors=True,
        )
        .filter(pl.col("feature") == "transcript")
        .select(
            region=pl.col("attributes").str.extract(r'gene_id "([^"]+)"')
                  .str.split(".").list.first(),
            enst_id=pl.col("attributes").str.extract(r'transcript_id "([^"]+)"')
                  .str.split(".").list.first(),
        )
        .filter(pl.col("region").is_not_null() & pl.col("enst_id").is_not_null())
        .unique()
    )


def step1b_revel(annos: pl.LazyFrame, revel_path: str, work_dir: str) -> pl.LazyFrame:
    """
    Add REVEL scores joined on (region, chrom, pos, ref, alt).
    REVEL is transcript-level (Ensembl_transcriptid column, semicolon-separated).
    We resolve transcript → region (ENSG) via a Gencode GTF (downloaded on
    demand) which covers every transcript including non-canonical ones —
    UniProt idmapping only covers UniProt-curated canonicals and was found to
    leave ~500k missense variants without a REVEL score.
    """
    logger.info("Step 1b: Merging REVEL scores")
    try:
        gtf_path = os.path.join(work_dir, "gencode.v40.annotation.gtf.gz")
        if not os.path.exists(gtf_path):
            if not aria2c_download(GENCODE_GTF_URL, gtf_path):
                logger.warning("  Gencode GTF download failed; skipping REVEL")
                return annos
        enst_to_ensg = _build_enst_to_ensg_from_gtf(gtf_path)
        logger.info(
            f"  ENST→ENSG map (Gencode): {len(enst_to_ensg):,} entries"
        )

        rev = (
            pl.scan_csv(
                revel_path,
                null_values=["."],
                schema_overrides={"chr": pl.Utf8},
                ignore_errors=True,
            )
            .with_columns(
                chrom="chr" + pl.col("chr").str.replace(r"^chr", "")
            )
            .rename({"grch38_pos": "pos", "REVEL": "revel_score"})
            .with_columns(
                pl.col("Ensembl_transcriptid").str.split(";")
            )
            .explode("Ensembl_transcriptid")
            .filter(
                pl.col("Ensembl_transcriptid").is_not_null()
                & (pl.col("Ensembl_transcriptid") != "")
            )
            .with_columns(
                enst_id=pl.col("Ensembl_transcriptid").str.split(".").list.first()
            )
            .select(["chrom", "pos", "ref", "alt", "revel_score", "enst_id"])
        )
        rev_with_region = (
            rev.join(enst_to_ensg.lazy(), on="enst_id", how="inner")
            .group_by(["region", "chrom", "pos", "ref", "alt"])
            .agg(pl.col("revel_score").max())
            .collect()  # small lookup table; collect before joining into main frame
        )
        annos = _idempotent_join(
            annos, rev_with_region.lazy(),
            on=["region", "chrom", "pos", "ref", "alt"],
        )
        logger.info("  REVEL OK")
    except Exception as e:
        logger.warning(f"  REVEL failed: {e}")
    return annos


def step1c_clinpred(annos: pl.LazyFrame, clinpred_path: str) -> pl.LazyFrame:
    """Add ClinPred pathogenicity scores (hg38 coordinates; no liftover needed)."""
    logger.info("Step 1c: Merging ClinPred scores")
    try:
        cp = (
            pl.scan_csv(
                clinpred_path,
                separator="\t",
                schema_overrides={"Chr": pl.Utf8, "ClinPred_Score": pl.Float32},
                ignore_errors=True,
            )
            .with_columns(
                chrom="chr" + pl.col("Chr").cast(pl.Utf8).str.replace(r"^chr", "")
            )
            .rename({"Start": "pos", "Ref": "ref", "Alt": "alt", "ClinPred_Score": "clinpred_score"})
            .select(["chrom", "pos", "ref", "alt", "clinpred_score"])
            .unique(subset=["chrom", "pos", "ref", "alt"], keep="first")
        )
        annos = _idempotent_join(annos, cp, on=["chrom", "pos", "ref", "alt"])
        logger.info("  ClinPred OK")
    except Exception as e:
        logger.warning(f"  ClinPred failed: {e}")
    return annos


def step1d_bayesdel(annos: pl.LazyFrame, bayesdel_path: str) -> pl.LazyFrame:
    """Add BayesDel noAF scores (joined on chrom, pos, ref, alt)."""
    logger.info("Step 1d: Merging BayesDel scores")
    try:
        bd = (
            pl.scan_csv(
                bayesdel_path,
                separator="\t",
                schema_overrides={"#Chr": pl.Utf8},
                ignore_errors=True,
            )
            .with_columns(
                chrom="chr" + pl.col("#Chr").cast(pl.Utf8).str.replace(r"^chr", "")
            )
            .rename({"Start": "pos", "BayesDel_nsfp33a_noAF": "bayes_del"})
            .select(["chrom", "pos", "ref", "alt", "bayes_del"])
            .unique(subset=["chrom", "pos", "ref", "alt"], keep="first")
        )
        annos = _idempotent_join(annos, bd, on=["chrom", "pos", "ref", "alt"])
        logger.info("  BayesDel OK")
    except Exception as e:
        logger.warning(f"  BayesDel failed: {e}")
    return annos


def step1e_popeve(annos: pl.LazyFrame, popeve_path: str) -> pl.LazyFrame:
    """Add popEVE / EVE / ESM1v scores from a VCF (joined on chrom, pos, ref, alt).

    Output columns: popeve, eve, esm1v, pop_adjusted_eve, pop_adjusted_esm1v
    """
    logger.info("Step 1e: Merging popEVE scores")
    try:
        vcf = (
            pl.scan_csv(
                popeve_path,
                separator="\t",
                comment_prefix="##",
                has_header=True,
                schema_overrides={
                    "#CHROM": pl.Utf8,
                    "POS": pl.Int64,
                    "REF": pl.Utf8,
                    "ALT": pl.Utf8,
                    "INFO": pl.Utf8,
                },
            )
            .select(
                chrom=pl.lit("chr")
                + pl.col("#CHROM").cast(pl.Utf8).str.replace(r"^chr", ""),
                pos=pl.col("POS"),
                ref=pl.col("REF"),
                alt=pl.col("ALT"),
                # Use (?:^|;) anchors so standalone EVE= / ESM1v= don't match the
                # pop-adjusted_ prefixed variants that also contain those substrings.
                popeve=pl.col("INFO")
                .str.extract(r"popEVE=([^;]+)")
                .cast(pl.Float32),
                eve=pl.col("INFO")
                .str.extract(r"(?:^|;)EVE=([^;]+)")
                .cast(pl.Float32),
                esm1v=pl.col("INFO")
                .str.extract(r"(?:^|;)ESM1v=([^;]+)")
                .cast(pl.Float32),
                pop_adjusted_eve=pl.col("INFO")
                .str.extract(r"pop-adjusted_EVE=([^;]+)")
                .cast(pl.Float32),
                pop_adjusted_esm1v=pl.col("INFO")
                .str.extract(r"pop-adjusted_ESM1v=([^;]+)")
                .cast(pl.Float32),
            )
            .unique(subset=["chrom", "pos", "ref", "alt"], keep="first")
        )
        annos = _idempotent_join(annos, vcf, on=["chrom", "pos", "ref", "alt"])
        logger.info("  popEVE OK")
    except Exception as e:
        logger.warning(f"  popEVE failed: {e}")
    return annos


def step1f_cpt1(
    annos: pl.LazyFrame,
    cpt1_dir: str,
    uniprot_map_df: pl.DataFrame,
) -> pl.LazyFrame:
    """
    Add CPT-1 per-protein scores (joined on region, amino_acids, prot_pos).
    cpt1_dir should contain a CPT1_all_proteins/ subdirectory with per-protein
    {entry_name}_HUMAN.csv.gz files downloaded from Zenodo record 8140323.
    """
    logger.info("Step 1f: Merging CPT-1 scores")
    schema_names = set(annos.collect_schema().names())
    if "protein_position" not in schema_names or "amino_acids" not in schema_names:
        logger.warning("  protein_position/amino_acids absent; skipping CPT-1")
        return annos
    if uniprot_map_df.is_empty():
        logger.warning("  UniProt map empty; skipping CPT-1")
        return annos
    try:
        # Build entry_name (gene symbol) → region mapping
        entry_to_region = (
            uniprot_map_df.select(["region", "entry_name"])
            .drop_nulls("entry_name")
            .unique()
        )
        present_regions = annos.select("region").unique().collect()
        cpt1_entry_to_region = (
            entry_to_region
            .join(present_regions, on="region", how="semi")
        )
        needed_entries = set(cpt1_entry_to_region["entry_name"].to_list())
        logger.info(f"  CPT-1: loading scores for {len(needed_entries)} proteins")

        # Build a lookup from UniProt Entry Name → file path, scanning all
        # subdirectories (CPT1_score_EVE_set, CPT1_score_no_EVE_set_*, etc.)
        cpt1_file_index = {
            fpath.name.replace(".csv.gz", ""): str(fpath)
            for fpath in Path(cpt1_dir).rglob("*_HUMAN.csv.gz")
        }
        logger.info(f"  CPT-1: index has {len(cpt1_file_index)} files in {cpt1_dir}")

        cpt1_dfs = []
        for entry_name in sorted(needed_entries):
            path = cpt1_file_index.get(entry_name)
            if path is None:
                continue
            try:
                df = (
                    pl.read_csv(path)
                    .rename({"CPT1_score": "cpt1_llr"})
                    .with_columns(
                        ref_aa=pl.col("mutant").str.slice(0, 1),
                        alt_aa=pl.col("mutant").str.slice(-1, 1),
                        prot_pos=pl.col("mutant").str.slice(
                            1, pl.col("mutant").str.len_chars() - 2
                        ),
                    )
                    .with_columns(
                        amino_acids=pl.col("ref_aa") + "/" + pl.col("alt_aa"),
                        entry_name=pl.lit(entry_name),
                    )
                    .select(["entry_name", "amino_acids", "prot_pos", "cpt1_llr"])
                )
                cpt1_dfs.append(df)
            except Exception as ef:
                logger.warning(f"  CPT-1: failed for {entry_name}: {ef}")

        if not cpt1_dfs:
            logger.warning("  CPT-1: no files loaded")
            return annos

        cpt1_scores = (
            pl.concat(cpt1_dfs)
            .join(cpt1_entry_to_region, on="entry_name", how="inner")
            .select(["region", "amino_acids", "prot_pos", "cpt1_llr"])
        )
        annos = (
            annos.with_columns(
                _prot_pos=pl.col("protein_position").str.split("/").list.get(0)
            )
            .join(
                cpt1_scores.lazy(),
                left_on=["region", "amino_acids", "_prot_pos"],
                right_on=["region", "amino_acids", "prot_pos"],
                how="left",
            )
            .drop("_prot_pos")
        )
        logger.info("  CPT-1 OK")
    except Exception as e:
        logger.warning(f"  CPT-1 failed: {e}")
    return annos



# ── Group 2: CADD / conservation / expression ──
def _tabix_fetch_scores(
    url: str,
    variants: pl.DataFrame,
    score_col_hint: str,
    out_col: str,
    chrom_prefix: str = "",
) -> pl.DataFrame:
    """
    Fetch scores for a set of variants from a remote tabix-indexed TSV/VCF.
    Groups variants by chromosome to minimise the number of tabix calls.
    Returns a DataFrame with columns: chrom, pos, ref, alt, {out_col}.

    chrom_prefix: prefix expected by the tabix index (e.g. "chr" for files that
    use "chr1", or "" for CADD-style "1"). Output chrom is always normalised
    back to "chr1" form regardless.
    """
    results = []
    # Resolve schema once from a tiny known region on chrom 1
    header: list[str] | None = None
    col_lower: list[str] | None = None
    probe_region = f"{chrom_prefix}1:925000-926000"
    try:
        proc0 = subprocess.run(
            [_tabix_bin(), "-h", url, probe_region],
            capture_output=True, text=True, check=True, timeout=60,
        )
        raw0 = proc0.stdout.splitlines()
        header_line = next((l for l in reversed(raw0) if l.startswith("#")), None)
        if header_line:
            header = header_line.lstrip("#").strip().split("\t")
            col_lower = [c.lower() for c in header]
    except Exception as e:
        logger.warning(f"  {out_col}: failed to fetch header: {e}")
        return pl.DataFrame(schema={"chrom": pl.Utf8, "pos": pl.Int64, "ref": pl.Utf8, "alt": pl.Utf8, out_col: pl.Float32})

    # Fallback for headerless tabix files (CADD, GPN-MSA): standard VCF-like format
    if col_lower is None:
        col_lower = ["chrom", "pos", "ref", "alt", "score"]

    chroms = variants["chrom"].unique().to_list()
    n_chroms = len(chroms)
    t0 = time.time()

    def _fetch_one(chrom_val: str, group: pl.DataFrame) -> pl.DataFrame | None:
        chrom_bare = str(chrom_val).lstrip("chr")
        chrom_query = f"{chrom_prefix}{chrom_bare}"
        positions = sorted(set(group["pos"].to_list()))

        # Write positions to a temp file for awk to use as a lookup set.
        # We stream the full chromosome through tabix and filter with awk at
        # C speed — this avoids random BGZF seeks (which are catastrophically
        # slow on network filesystems for millions of query positions).
        pos_path = os.path.join(
            WORK_DIR, f"_tabix_pos_{chrom_bare}_{os.getpid()}.txt"
        )
        with open(pos_path, "w") as f:
            for p in positions:
                f.write(f"{p}\n")

        logger.info(
            f"  {out_col}: {chrom_val} — {len(positions):,} positions, streaming"
        )
        t_chrom = time.time()
        try:
            # tabix streams the full chromosome; awk keeps only matching positions
            tabix_proc = subprocess.Popen(
                [_tabix_bin(), url, f"{chrom_query}:1-999999999"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            awk_script = (
                f'BEGIN{{while((getline<"{pos_path}")>0)p[$0]=1}}'
                f'($2 in p)'
            )
            awk_proc = subprocess.Popen(
                ["awk", "-F", "\t", awk_script],
                stdin=tabix_proc.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            tabix_proc.stdout.close()  # let tabix get SIGPIPE if awk exits early
            out_bytes, _ = awk_proc.communicate(timeout=600)
            tabix_proc.wait(timeout=30)

            if tabix_proc.returncode not in (0, None):
                err = tabix_proc.stderr.read().decode(errors="replace").strip()
                logger.warning(f"  tabix stream failed for {chrom_val}: {err[:300]}")
                return None

            t_io = time.time() - t_chrom
            if not out_bytes:
                logger.info(f"  {out_col}: {chrom_val} — 0 hits ({t_io:.0f}s I/O)")
                return None
            # Use Polars CSV reader (C-speed) instead of Python decode+splitlines+split
            import io as _io
            df = pl.read_csv(
                _io.BytesIO(out_bytes),
                separator="\t",
                has_header=False,
                new_columns=col_lower,
                infer_schema_length=0,
            )
            if df.is_empty():
                logger.info(f"  {out_col}: {chrom_val} — 0 hits ({t_io:.0f}s I/O)")
                return None
            for alias, targets in [
                ("chrom", ("chr", "#chrom", "chromosome")),
                ("ref", ("#ref",)),
                ("alt", ("#alt",)),
            ]:
                for t in targets:
                    if t in df.columns and alias not in df.columns:
                        df = df.rename({t: alias})
                        break
            df = df.with_columns(
                pl.when(pl.col("chrom").str.starts_with("chr"))
                .then(pl.col("chrom"))
                .otherwise(pl.lit("chr") + pl.col("chrom"))
                .alias("chrom")
            )
            key_cols = {"chrom", "pos", "ref", "alt"}
            score_col = next(
                (c for c in df.columns if score_col_hint in c),
                next((c for c in reversed(df.columns) if c not in key_cols), None),
            )
            if score_col is None:
                return None
            df = (
                df.with_columns(
                    pl.col("pos").cast(pl.Int64),
                    pl.col(score_col).cast(pl.Float32).alias(out_col),
                )
                .select(["chrom", "pos", "ref", "alt", out_col])
                .join(group.select(["chrom", "pos", "ref", "alt"]), on=["chrom", "pos", "ref", "alt"], how="semi")
            )
            t_total = time.time() - t_chrom
            logger.info(
                f"  {out_col}: {chrom_val} done ({len(df):,} hits, {t_io:.0f}s I/O + {t_total-t_io:.0f}s parse, {t_total:.0f}s total)"
            )
            return df
        except Exception as e:
            logger.warning(f"  tabix fetch failed for {chrom_val}: {e}")
            return None
        finally:
            try:
                os.unlink(pos_path)
            except OSError:
                pass

    groups = [(str(cv), g) for (cv,), g in variants.group_by(["chrom"])]
    logger.info(f"  {out_col}: dispatching {len(groups)} chromosomes in parallel ...")
    with ThreadPoolExecutor(max_workers=min(8, len(groups))) as ex:
        futs = {ex.submit(_fetch_one, cv, g): cv for cv, g in groups}
        completed = 0
        for fut in as_completed(futs):
            completed += 1
            df = fut.result()
            if df is not None:
                results.append(df)
            logger.info(
                f"  {out_col}: progress {completed}/{n_chroms} chroms, "
                f"{time.time()-t0:.0f}s total"
            )

    if not results:
        return pl.DataFrame(schema={"chrom": pl.Utf8, "pos": pl.Int64, "ref": pl.Utf8, "alt": pl.Utf8, out_col: pl.Float32})
    return pl.concat(results, how="diagonal")


def step2a_cadd(annos: pl.LazyFrame) -> pl.LazyFrame:
    """
    Add CADD raw scores via remote tabix queries against the CADD GRCh38 server.
    SNVs use whole_genome_SNVs.tsv.gz; indels use the gnomAD v4.0 genomes indel
    file (CADD publishes no whole-genome indel score set, only this
    gnomAD-variant one — so indel coverage is limited to variants gnomAD has
    scored, unlike the SNV file's full-genome coverage).
    Output column: cadd_raw.
    """
    logger.info("Step 2a: Merging CADD scores (local tabix)")
    try:
        variants = annos.select(["chrom", "pos", "ref", "alt"]).unique().collect()

        # Split SNVs and indels
        snvs = variants.filter(
            (pl.col("ref").str.len_bytes() == 1) & (pl.col("alt").str.len_bytes() == 1)
        )
        indels = variants.filter(
            (pl.col("ref").str.len_bytes() != 1) | (pl.col("alt").str.len_bytes() != 1)
        )

        frames = []
        if not snvs.is_empty():
            cadd_local = _ensure_local_tabix(CADD_SNV_URL)
            if cadd_local is None:
                logger.warning("  CADD: SNV download failed, skipping SNVs")
            else:
                frames.append(_tabix_fetch_scores(cadd_local, snvs, "rawscore", "cadd_raw"))

        if not indels.is_empty():
            logger.info(f"  CADD: {len(indels):,} indels — querying gnomAD indel file")
            cadd_indel_local = _ensure_local_tabix(CADD_INDEL_URL)
            if cadd_indel_local is None:
                logger.warning(f"  CADD: indel download failed, skipping {len(indels):,} indels")
            else:
                frames.append(_tabix_fetch_scores(cadd_indel_local, indels, "rawscore", "cadd_raw"))

        if not frames:
            logger.warning("  CADD: no scores retrieved")
            return annos

        cadd_df = pl.concat(frames, how="diagonal")
        annos = _idempotent_join(annos, cadd_df.lazy(), on=["chrom", "pos", "ref", "alt"])
        logger.info("  CADD OK")
    except Exception as e:
        logger.warning(f"  CADD failed: {e}")
    return annos


def step2b_gpn_msa(annos: pl.LazyFrame) -> pl.LazyFrame:
    """
    Add GPN-MSA conservation scores via remote tabix against the HuggingFace bgz.
    No full download required; only variant-overlapping rows are fetched.
    Requires tabix on PATH. Output column: gpn_score.
    """
    logger.info("Step 2b: Merging GPN-MSA scores (local tabix)")
    try:
        variants = annos.select(["chrom", "pos", "ref", "alt"]).unique().collect()
        gpn_local = _ensure_local_tabix(GPN_MSA_URL)
        if gpn_local is None:
            logger.warning("  GPN-MSA: download failed, skipping")
            return annos
        gpn_df = _tabix_fetch_scores(
            gpn_local, variants, "score", "gpn_score", chrom_prefix=""
        )
        if gpn_df.is_empty():
            logger.warning("  GPN-MSA: no scores retrieved")
            return annos
        annos = _idempotent_join(annos, gpn_df.lazy(), on=["chrom", "pos", "ref", "alt"])
        logger.info("  GPN-MSA OK")
    except Exception as e:
        logger.warning(f"  GPN-MSA failed: {e}")
    return annos


def _bigwig_fetch_scores(
    url: str,
    variants: pl.DataFrame,
    out_col: str,
    local_path: str | None = None,
) -> pl.DataFrame:
    """
    Fetch per-position scores from a bigWig file using pyBigWig.
    Downloads the file locally first (via aria2c), then queries in bulk:
    one bw.values() call per chromosome covering the full variant range,
    with numpy indexing to extract exact positions — no per-variant HTTP calls.
    local_path: where to cache the downloaded file; defaults to WORK_DIR/<filename>.
    """
    try:
        import pyBigWig  # noqa: PLC0415
    except ImportError:
        logger.warning(
            f"  pyBigWig not installed; skipping {out_col}. "
            "Add pyBigWig to your conda environment."
        )
        return pl.DataFrame(schema={"chrom": pl.Utf8, "pos": pl.Int64, "ref": pl.Utf8, "alt": pl.Utf8, out_col: pl.Float32})

    empty = pl.DataFrame(schema={"chrom": pl.Utf8, "pos": pl.Int64, "ref": pl.Utf8, "alt": pl.Utf8, out_col: pl.Float32})

    # Download locally so we can do bulk per-chromosome queries without
    # making one HTTP call per variant.
    if local_path is None:
        local_path = os.path.join(WORK_DIR, os.path.basename(url.split("?")[0]))
    if not os.path.exists(local_path):
        logger.info(f"  Downloading bigWig: {url}")
        ok = aria2c_download(url, local_path)
        if not ok:
            logger.warning(f"  bigWig download failed; skipping {out_col}")
            return empty

    frames = []
    n_chroms = variants["chrom"].n_unique()
    t0 = time.time()
    # Cap peak per-window allocation: 25 Mb × float32 ≈ 100 MB/window.
    WINDOW = 25_000_000
    try:
        bw = pyBigWig.open(local_path)
        for done, ((chrom_val,), group) in enumerate(variants.group_by(["chrom"]), 1):
            chrom = str(chrom_val)
            positions = group["pos"].to_numpy()  # 1-based VCF positions
            min_pos = int(positions.min())
            max_pos = int(positions.max())
            scores = np.full(len(positions), np.nan, dtype=np.float32)
            n_windows = ((max_pos - min_pos) // WINDOW) + 1
            try:
                for wi, win_start in enumerate(
                    range(min_pos, max_pos + 1, WINDOW), 1
                ):
                    win_end = min(win_start + WINDOW, max_pos + 1)
                    mask = (positions >= win_start) & (positions < win_end)
                    if not mask.any():
                        continue
                    # bw.values uses 0-based half-open [start, end)
                    vals = bw.values(chrom, win_start - 1, win_end - 1, numpy=True)
                    offsets = positions[mask] - win_start
                    scores[mask] = vals[offsets]
                    logger.info(
                        f"  {out_col}: {chrom} window {wi}/{n_windows} "
                        f"({mask.sum():,} positions, {time.time()-t0:.0f}s elapsed)"
                    )
            except Exception as we:
                logger.warning(f"  bigWig window error for {chrom}: {we}")
            frames.append(pl.DataFrame({
                "chrom": [chrom] * len(positions),
                "pos": positions.tolist(),
                "ref": group["ref"].to_list(),
                "alt": group["alt"].to_list(),
                out_col: scores.tolist(),
            }))
            elapsed = time.time() - t0
            rate = done / elapsed if elapsed > 0 else 0
            eta = (n_chroms - done) / rate if rate > 0 else float("inf")
            logger.info(
                f"  {out_col}: {chrom} done ({done}/{n_chroms} chroms, "
                f"{elapsed:.0f}s elapsed, ETA {eta:.0f}s)"
            )
        bw.close()
    except Exception as e:
        logger.warning(f"  bigWig fetch error for {out_col}: {e}")
        return empty

    if not frames:
        return empty
    return pl.concat(frames, how="vertical").with_columns(
        pl.when(
            (pl.col("ref").str.len_bytes() == 1) &
            (pl.col("alt").str.len_bytes() == 1) &
            ~pl.col(out_col).is_nan()
        )
        .then(pl.col(out_col).cast(pl.Float32))
        .otherwise(pl.lit(None, dtype=pl.Float32))
        .alias(out_col)
    )


def step2c1_phylop(annos: pl.LazyFrame) -> pl.LazyFrame:
    """
    Add PhyloP 100-way vertebrate conservation scores.
    Downloads the bigWig locally then queries per chromosome in bulk.
    Output column: phylop_100way.
    """
    logger.info("Step 2c1: Merging PhyloP 100-way scores (bigWig, local download)")
    try:
        variants = annos.select(["chrom", "pos", "ref", "alt"]).unique().collect()
        scores_df = _bigwig_fetch_scores(
            PHYLOP_100WAY_URL, variants, "verphylop",
            local_path=f"{WORK_DIR}/hg38.phyloP100way.bw",
        )
        if scores_df.is_empty():
            logger.warning("  PhyloP 100-way: no scores retrieved (pyBigWig unavailable or download failed)")
        else:
            annos = _idempotent_join(annos, scores_df.lazy(), on=["chrom", "pos", "ref", "alt"])
            logger.info("  PhyloP 100-way OK")
    except Exception as e:
        logger.warning(f"  PhyloP 100-way failed: {e}")
    return annos


def step2c2_phylop_mammalian(annos: pl.LazyFrame) -> pl.LazyFrame:
    """Add PhyloP 17-way mammalian conservation scores (mamphylop)."""
    logger.info("Step 2c2: Merging PhyloP 17-way mammalian scores (bigWig, local download)")
    try:
        variants = annos.select(["chrom", "pos", "ref", "alt"]).unique().collect()
        scores_df = _bigwig_fetch_scores(
            PHYLOP_17WAY_URL, variants, "mamphylop",
            local_path=f"{WORK_DIR}/hg38.phyloP17way.bw",
        )
        if scores_df.is_empty():
            logger.warning("  PhyloP 17-way: no scores retrieved (pyBigWig unavailable or download failed)")
        else:
            annos = _idempotent_join(annos, scores_df.lazy(), on=["chrom", "pos", "ref", "alt"])
            logger.info("  PhyloP 17-way mammalian OK")
    except Exception as e:
        logger.warning(f"  PhyloP 17-way mammalian failed: {e}")
    return annos


def step2c3_phylop_primate(annos: pl.LazyFrame) -> pl.LazyFrame:
    """Add PhyloP 30-way primate conservation scores (priphylop)."""
    logger.info("Step 2c3: Merging PhyloP 30-way primate scores (bigWig, local download)")
    try:
        variants = annos.select(["chrom", "pos", "ref", "alt"]).unique().collect()
        scores_df = _bigwig_fetch_scores(
            PHYLOP_30WAY_URL, variants, "priphylop",
            local_path=f"{WORK_DIR}/hg38.phyloP30way.bw",
        )
        if scores_df.is_empty():
            logger.warning("  PhyloP 30-way: no scores retrieved (pyBigWig unavailable or download failed)")
        else:
            annos = _idempotent_join(annos, scores_df.lazy(), on=["chrom", "pos", "ref", "alt"])
            logger.info("  PhyloP 30-way primate OK")
    except Exception as e:
        logger.warning(f"  PhyloP 30-way primate failed: {e}")
    return annos



# ── Group 3: Derived columns + fill nulls ──
def step3a_derived_columns(annos: pl.LazyFrame) -> pl.LazyFrame:
    """Add columns derived from already-merged annotations."""
    logger.info("Step 3a: Adding derived columns")
    schema = set(annos.collect_schema().names())
    exprs = []

    # ENCODE derived columns
    encode_mapping = {
        ("encode_pls",): "core_promoter",
        ("encode_pels",): "proximal_promoter",
        ("encode_dels",): "encode_enhancer",
    }
    for (src,), dst in encode_mapping.items():
        if src in schema:
            exprs.append(pl.col(src).cast(pl.Int8).alias(dst))
    if "encode_pls" in schema and "encode_pels" in schema:
        exprs.append(
            (pl.col("encode_pls") | pl.col("encode_pels"))
            .cast(pl.Int8)
            .alias("encode_promoter")
        )
    if "encode_tf" in schema and "encode_ca_tf" in schema:
        exprs.append(
            (pl.col("encode_tf") | pl.col("encode_ca_tf"))
            .cast(pl.Int8)
            .alias("encode_any_tf")
        )

    # Indel length categories
    if "variant_length" in schema and "is_insertion" in schema and "is_deletion" in schema:
        exprs.extend(
            [
                ((pl.col("is_deletion") == 1) & (pl.col("variant_length") == 2))
                .cast(pl.Int8)
                .alias("1bp_del"),
                ((pl.col("is_insertion") == 1) & (pl.col("variant_length") == 2))
                .cast(pl.Int8)
                .alias("1bp_ins"),
                (
                    (pl.col("is_deletion") == 1)
                    & pl.col("variant_length").is_between(3, 6)
                )
                .cast(pl.Int8)
                .alias("2_5bp_del"),
                (
                    (pl.col("is_insertion") == 1)
                    & pl.col("variant_length").is_between(3, 6)
                )
                .cast(pl.Int8)
                .alias("2_5bp_ins"),
                ((pl.col("is_deletion") == 1) & (pl.col("variant_length") > 6))
                .cast(pl.Int8)
                .alias("gt_5bp_del"),
                ((pl.col("is_insertion") == 1) & (pl.col("variant_length") > 6))
                .cast(pl.Int8)
                .alias("gt_5bp_ins"),
            ]
        )

    if exprs:
        annos = annos.with_columns(exprs)
    return annos


def step3b_fill_nulls(
    annos: pl.LazyFrame, fill_defaults: dict
) -> pl.LazyFrame:
    """
    For each annotation column with a fill default:
      1. Add {col}_is_na (Int8: 1 if null before filling)
      2. Fill null with the default value
    Consequence columns (consequence_*, five_prime_utr_*) are filled with 0,
    no _is_na indicator.
    Cast all float columns to Float32.
    """
    logger.info("Step 3b: Filling nulls and adding _is_na indicators")
    schema = annos.collect_schema()

    # Consequence columns → fill with 0, no _is_na
    consequence_cols = [
        c
        for c in schema.names()
        if c.startswith("consequence_") or c.startswith("five_prime_utr_")
    ]
    if consequence_cols:
        annos = annos.with_columns(
            [pl.col(c).fill_null(0) for c in consequence_cols]
        )

    # Cast float64 → float32
    float64_cols = [c for c, t in schema.items() if t == pl.Float64]
    if float64_cols:
        annos = annos.with_columns(
            [pl.col(c).cast(pl.Float32) for c in float64_cols]
        )

    # For each column with a fill default: add _is_na, then fill
    schema_after = annos.collect_schema()
    schema_names = set(schema_after.names())
    for col, fill_val in fill_defaults.items():
        if col in schema_names:
            col_expr = pl.col(col)
            # Boolean columns (e.g. encode_* cCRE flags) can't be filled with an
            # int default directly ([bool, int] is ambiguous). Cast to Int8 first
            # so null→0 and True/False→1/0, matching step4c_encode's encoding.
            if schema_after[col] == pl.Boolean and isinstance(fill_val, int):
                col_expr = col_expr.cast(pl.Int8)
            annos = annos.with_columns(
                pl.col(col).is_null().cast(pl.Int8).alias(f"{col}_is_na")
            ).with_columns(col_expr.fill_null(fill_val))

    return annos

# ── Group 4: Post-fill protein-domain & external annotations ──
def step4a_plddt(
    annos: pl.LazyFrame, uniprot_map: pl.DataFrame, work_dir: str
) -> pl.LazyFrame:
    """Add AlphaFold pLDDT per-residue confidence (Float32, 0–100).

    Streams per-protein PDB.gz members from the AlphaFold human-proteome tar
    (only the proteins present in `annos`), parses the B-factor column on CA
    atoms, then joins on (region, prot_pos) via uniprot_map.
    """
    logger.info("Step 4a: Adding AlphaFold pLDDT scores")

    schema_names = set(annos.collect_schema().names())
    if "protein_position" not in schema_names or "region" not in schema_names:
        logger.warning("  protein_position/region absent; skipping pLDDT")
        return annos

    present_regions = annos.select("region").unique().collect()
    region_to_uniprot = (
        uniprot_map.select(["region", "uniprot_id"])
        .drop_nulls("uniprot_id")
        .join(present_regions, on="region", how="semi")
    )
    needed_uniprots = set(region_to_uniprot["uniprot_id"].to_list())
    logger.info(f"  Need pLDDT for {len(needed_uniprots):,} proteins")

    tar_path = os.path.join(work_dir, "UP000005640_9606_HUMAN_v6.tar")
    if not os.path.exists(tar_path):
        if not aria2c_download(ALPHAFOLD_HUMAN_TAR_URL, tar_path):
            logger.warning("  AlphaFold tar download failed; skipping pLDDT")
            return annos

    plddt_records = []
    with tarfile.open(tar_path, "r") as tf:
        for member in tf:
            if not member.isfile():
                continue
            name = member.name
            if not name.startswith("AF-") or not name.endswith(".pdb.gz"):
                continue
            uid = name.split("-")[1]
            if uid not in needed_uniprots:
                continue
            f = tf.extractfile(member)
            if f is None:
                continue
            content = gzip.decompress(f.read()).decode("ascii", errors="ignore")
            for line in content.splitlines():
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    try:
                        res_idx = int(line[22:26].strip())
                        plddt = float(line[60:66].strip())
                    except ValueError:
                        continue
                    plddt_records.append((uid, res_idx, plddt))

    if not plddt_records:
        logger.warning("  No pLDDT records extracted; skipping")
        return annos

    plddt_df = pl.DataFrame(
        plddt_records,
        schema={"uniprot_id": pl.Utf8, "residue_number": pl.Int32, "plddt": pl.Float32},
        orient="row",
    )
    logger.info(
        f"  Extracted pLDDT for {plddt_df['uniprot_id'].n_unique():,} proteins, "
        f"{plddt_df.shape[0]:,} residues"
    )

    plddt_with_region = (
        plddt_df.join(region_to_uniprot, on="uniprot_id", how="inner")
        .group_by(["region", "residue_number"])
        .agg(pl.col("plddt").mean())
    )

    annos = (
        annos.with_columns(
            prot_pos=pl.col("protein_position").str.extract(r"(\d+)", 0).cast(pl.Int32)
        )
        .join(
            plddt_with_region.lazy(),
            left_on=["region", "prot_pos"],
            right_on=["region", "residue_number"],
            how="left",
        )
        .drop("prot_pos")
    )
    logger.info("  pLDDT OK")
    return annos


def step4b_pioneer_interface(
    annos: pl.LazyFrame, uniprot_map: pl.DataFrame, work_dir: str
) -> pl.LazyFrame:
    """Add PIONEER High-confidence interface residue flag (Int8, 0/1).

    Downloads the PIONEER human predictions file (~173 MB), explodes the
    bracketed residue lists for each interaction, and left-joins on
    (region, prot_pos) via uniprot_map.  Nulls are filled with 0.
    """
    logger.info("Step 4b: Adding PIONEER interface residues")

    schema_names = set(annos.collect_schema().names())
    if "protein_position" not in schema_names or "region" not in schema_names:
        logger.warning(
            "  protein_position/region absent; setting is_pioneer_interface_high=0"
        )
        return annos.with_columns(
            is_pioneer_interface_high=pl.lit(0, dtype=pl.Int8)
        )

    pioneer_path = os.path.join(work_dir, "pioneer_high_human.txt")
    if not aria2c_download(PIONEER_HIGH_URL, pioneer_path):
        logger.warning(
            "  PIONEER download failed; setting is_pioneer_interface_high=0"
        )
        return annos.with_columns(
            is_pioneer_interface_high=pl.lit(0, dtype=pl.Int8)
        )

    present_regions = annos.select("region").unique().collect()
    region_to_uniprot = (
        uniprot_map.select(["region", "uniprot_id"])
        .drop_nulls("uniprot_id")
        .join(present_regions, on="region", how="semi")
    )
    needed_uniprots = set(region_to_uniprot["uniprot_id"].to_list())

    def parse_residues(s: str) -> list[int]:
        s = s.strip("[]")
        if not s:
            return []
        out = []
        for x in s.split(","):
            x = x.strip()
            if x:
                try:
                    out.append(int(x))
                except ValueError:
                    pass
        return out

    records: list[tuple[str, int]] = []
    with open(pioneer_path) as f:
        next(f, None)  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            uid1, uid2, res1_str, res2_str = parts[0], parts[1], parts[3], parts[4]
            if uid1 in needed_uniprots:
                for pos in parse_residues(res1_str):
                    records.append((uid1, pos))
            if uid2 in needed_uniprots and uid2 != uid1:
                for pos in parse_residues(res2_str):
                    records.append((uid2, pos))

    if not records:
        logger.warning("  No PIONEER records matched; setting is_pioneer_interface_high=0")
        return annos.with_columns(
            is_pioneer_interface_high=pl.lit(0, dtype=pl.Int8)
        )

    pioneer_df = (
        pl.DataFrame(
            records,
            schema={"uniprot_id": pl.Utf8, "residue_number": pl.Int32},
            orient="row",
        )
        .unique()
        .join(region_to_uniprot, on="uniprot_id", how="inner")
        .select(["region", "residue_number"])
        .unique()
        .with_columns(is_pioneer_interface_high=pl.lit(1, dtype=pl.Int8))
    )
    logger.info(
        f"  PIONEER interface residues: {pioneer_df.shape[0]:,} pairs across "
        f"{pioneer_df['region'].n_unique():,} genes"
    )

    annos = (
        annos.with_columns(
            prot_pos=pl.col("protein_position").str.extract(r"(\d+)", 0).cast(pl.Int32)
        )
        .join(
            pioneer_df.lazy(),
            left_on=["region", "prot_pos"],
            right_on=["region", "residue_number"],
            how="left",
        )
        .with_columns(
            is_pioneer_interface_high=pl.col("is_pioneer_interface_high").fill_null(0)
        )
        .drop("prot_pos")
    )
    logger.info("  PIONEER OK")
    return annos


def step4c_encode(annos: pl.LazyFrame, encode_bed_path: str) -> pl.LazyFrame:
    """Add ENCODE cCRE annotations via per-chromosome interval overlap (join_where).

    cCREs of different types can overlap, so a single asof join is not sufficient.
    join_where (polars IEJOIN) correctly finds all overlapping intervals per variant.
    Processed per chromosome to keep peak memory reasonable.
    """
    logger.info("Step 4c: Merging ENCODE cCRE annotations")
    try:
        encode_bool_cols = list(CCRE_COLUMN_MAP.values())

        # WENG Lab Registry V4: chrom, start (0-based BED), end, elem_acc, cre_acc, feature
        # feature contains compound values like "pELS,CTCF-bound" — extract primary type
        cres = pl.read_csv(
            encode_bed_path,
            separator="\t",
            has_header=False,
            new_columns=["chrom", "start", "end", "elem_acc", "cre_acc", "feature"],
            schema_overrides={"chrom": pl.Utf8, "start": pl.Int64, "end": pl.Int64, "feature": pl.Utf8},
        ).select(["chrom", "start", "end", "feature"]).with_columns(
            pl.col("feature").str.split(",").list.first().alias("feature")
        )

        # Convert VCF 1-based pos to 0-based for BED comparison
        variants_pos = (
            annos.select(["id", "chrom", "pos"]).unique().collect()
            .with_columns((pl.col("pos") - 1).alias("pos0"))
        )

        # Per-chromosome join_where: finds ALL overlapping cCREs (handles type overlaps)
        overlap_parts: list[pl.DataFrame] = []
        for chrom in variants_pos["chrom"].unique().to_list():
            v = variants_pos.filter(pl.col("chrom") == chrom).select(["id", "pos0"])
            c = cres.filter(pl.col("chrom") == chrom).select(["start", "end", "feature"])
            if v.is_empty() or c.is_empty():
                continue
            hits = (
                v.join_where(c, pl.col("pos0") >= pl.col("start"), pl.col("pos0") < pl.col("end"))
                .select(["id", "feature"])
            )
            if not hits.is_empty():
                overlap_parts.append(hits)

        overlaps = pl.concat(overlap_parts, how="vertical") if overlap_parts else pl.DataFrame({"id": pl.Series([], dtype=pl.Utf8), "feature": pl.Series([], dtype=pl.Utf8)})

        # One boolean flag per cCRE type; max aggregation handles multiple hits per variant
        feature_flags = (
            overlaps
            .with_columns([
                (pl.col("feature") == feat).cast(pl.Int8).alias(col)
                for feat, col in CCRE_COLUMN_MAP.items()
            ])
            .group_by("id")
            .agg([pl.col(c).max() for c in encode_bool_cols])
        )

        encode_pl = (
            variants_pos.select("id")
            .join(feature_flags, on="id", how="left")
            .with_columns([pl.col(c).fill_null(0) for c in encode_bool_cols])
            # V3 categories absent in V4 — add as zeros for downstream compatibility
            .with_columns([pl.lit(0, dtype=pl.Int8).alias(c) for c in CCRE_ZERO_COLS_V4])
        )

        annos = _idempotent_join(annos, encode_pl.lazy(), on=["id"], validate="m:1")
        logger.info("  ENCODE OK")
    except Exception as e:
        logger.warning(f"  ENCODE failed: {e}")
    return annos


def step4d_clinvar(annos: pl.LazyFrame, clinvar_vcf_path: str) -> pl.LazyFrame:
    """Add ClinVar clinical_significance and pathogenicity dummy columns."""
    logger.info("Step 4d: Merging ClinVar annotations")
    try:
        needed_ids = annos.select("id").unique().collect()
        clinvar = (
            pl.scan_csv(
                clinvar_vcf_path,
                separator="\t",
                comment_prefix="##",
                schema_overrides={"#CHROM": pl.Utf8, "POS": pl.Int64},
                ignore_errors=True,
            )
            .rename({"#CHROM": "chrom", "POS": "pos", "REF": "ref", "ALT": "alt"})
            .with_columns(
                chrom="chr"
                + pl.col("chrom").cast(pl.Utf8).str.replace(r"^chr", ""),
            )
            .with_columns(
                id=pl.col("chrom")
                + ":"
                + pl.col("pos").cast(pl.Utf8)
                + ":"
                + pl.col("ref")
                + ":"
                + pl.col("alt"),
                clinical_significance=pl.col("INFO").str.extract(
                    r"CLNSIG=([^;]+)", 1
                ),
            )
            .select(["id", "clinical_significance"])
            .join(needed_ids.lazy(), on="id", how="semi")
            .unique()
            .collect()
        )
        annos = _idempotent_join(annos, clinvar.lazy(), on=["id"], validate="m:1")
        # Pathogenicity dummies (str.contains handles multi-value CLNSIG strings)
        annos = annos.with_columns(
            clinvar_patho=pl.col("clinical_significance")
            .str.contains("^Pathogenic$|^Pathogenic/")
            .cast(pl.Int8)
            .fill_null(0),
            clinvar_likely_patho=pl.col("clinical_significance")
            .str.contains("Likely_pathogenic")
            .cast(pl.Int8)
            .fill_null(0),
            clinvar_benign=pl.col("clinical_significance")
            .str.contains("^Benign$|^Benign/")
            .cast(pl.Int8)
            .fill_null(0),
            clinvar_likely_benign=pl.col("clinical_significance")
            .str.contains("Likely_benign")
            .cast(pl.Int8)
            .fill_null(0),
        )
        logger.info("  ClinVar OK")
    except Exception as e:
        logger.warning(f"  ClinVar failed: {e}")
    return annos


# ── Group 0: VEP / LOFTEE post-processing (AoU-native columns) ──────────────
#
# Ported from the vep_loftee_parallel applet's post_process_vep() and
# add_vep_structural_features(), adapted to the AoU variant schema and rewritten
# polars-only (the applet used pandas + pyranges).
#
# Differences vs. the applet, driven by the AoU schema:
#   * LOFTEE fields (lof/lof_filter/lof_flags) are already parsed columns in the
#     AoU file, so the 'extra'-column extraction step is dropped.
#   * AoU has no CDS_position column. relative_cds_position is instead derived
#     from protein_start and the per-transcript CDS length read from the Gencode
#     GTF (an approximation: protein_start*3 / CDS_length_nt).
#   * next_in_frame_relative is omitted: it needs cds_position, strand and allele
#     (none present in AoU) plus a reference FASTA + pyranges.


def _load_gtf_polars(gtf_path: str) -> pl.DataFrame:
    """Read a (possibly gzipped) Gencode GTF into a polars DataFrame."""
    return pl.read_csv(
        gtf_path,
        separator="\t",
        comment_prefix="#",
        has_header=False,
        new_columns=[
            "chrom", "source", "feature", "start", "end",
            "score", "strand", "frame", "attributes",
        ],
        schema_overrides={
            "chrom": pl.Utf8, "start": pl.Int64, "end": pl.Int64,
            "feature": pl.Utf8, "strand": pl.Utf8, "attributes": pl.Utf8,
        },
        ignore_errors=True,
    )


def _gtf_gene_features(gtf: pl.DataFrame) -> pl.DataFrame:
    """Per-gene TSS, length and name from GTF 'gene' rows (protein-coding only).

    Returns columns: region, gene_name, gene_length, _tss, _gene_strand.
    """
    return (
        gtf.filter(pl.col("feature") == "gene")
        .with_columns(
            region=pl.col("attributes").str.extract(r'gene_id "([^"]+)"')
                  .str.split(".").list.first(),
            gene_name=pl.col("attributes").str.extract(r'gene_name "([^"]+)"'),
            gene_type=pl.col("attributes")
                  .str.extract(r'gene_(?:type|biotype) "([^"]+)"'),
        )
        .filter(pl.col("gene_type") == "protein_coding")
        .with_columns(
            gene_length=pl.col("end") - pl.col("start") + 1,
            _tss=pl.when(pl.col("strand") == "+")
                  .then(pl.col("start"))
                  .otherwise(pl.col("end")),
            _gene_strand=pl.col("strand"),
        )
        .select(["region", "gene_name", "gene_length", "_tss", "_gene_strand"])
        .unique(subset=["region"])
    )


def step0a_vep_loftee_dummies(
    annos: pl.LazyFrame, consequence_col: str = "consequence_terms"
) -> pl.LazyFrame:
    """LOFTEE one-hot dummies + VEP consequence one-hot encoding.

    LOFTEE: from the `lof` column → loftee_hc, loftee_lc (Int8), plus
    loftee_hc_is_na / loftee_lc_is_na.

    Consequence: from `consequence_col` (default 'consequence_terms', the
    per-transcript array extracted from Hail; falls back to 'vep_consequence')
    → one Int8 column per term, named consequence_{term}. Handles both native
    Hail list columns and VEP '&'-joined string columns. These columns are later
    filled with 0 by step3b_fill_nulls.
    """
    logger.info("Step 0a: LOFTEE dummies + consequence one-hot")
    schema = set(annos.collect_schema().names())

    # ── LOFTEE dummies ────────────────────────────────────────────────────
    if "lof" in schema:
        annos = annos.with_columns(
            loftee_hc=pl.col("lof").eq("HC").cast(pl.Int8).fill_null(0),
            loftee_lc=pl.col("lof").eq("LC").cast(pl.Int8).fill_null(0),
            loftee_hc_is_na=pl.col("lof").is_null().cast(pl.Int8),
            loftee_lc_is_na=pl.col("lof").is_null().cast(pl.Int8),
        )
        logger.info("  LOFTEE dummies OK")
    else:
        logger.warning("  'lof' column absent; skipping LOFTEE dummies")

    # ── Consequence one-hot ───────────────────────────────────────────────
    col = (
        consequence_col if consequence_col in schema
        else ("vep_consequence" if "vep_consequence" in schema else None)
    )
    if col is None:
        logger.warning("  No consequence column found; skipping one-hot")
        return annos

    # Strip array punctuation/quotes/spaces, normalise '&' → ',', then split.
    annos = annos.with_columns(
        _cons=pl.col(col).cast(pl.Utf8)
              .str.replace_all("[\"'\\[\\] ]", "")
              .str.replace_all("&", ",")
              .str.split(",")
    )
    terms = (
        annos.select(pl.col("_cons"))
        .explode("_cons")
        .drop_nulls("_cons")
        .filter(pl.col("_cons") != "")
        .unique()
        .collect()
        .get_column("_cons")
        .to_list()
    )
    if not terms:
        logger.warning("  No consequence terms parsed; skipping one-hot")
        return annos.drop("_cons")
    exprs = [
        pl.col("_cons").list.contains(t).cast(pl.Int8).alias(f"consequence_{t}")
        for t in sorted(terms)
    ]
    annos = annos.with_columns(exprs).drop("_cons")
    logger.info(f"  Consequence one-hot OK ({len(terms)} terms from '{col}')")
    return annos


def _next_inframe_atg_distance(seq: str, search_init: int = 3) -> int:
    """Return nt distance from position `search_init` to the next in-frame ATG.

    Searches codons starting at `search_init` (0-based, must be frame-0
    relative to the CDS start). Returns the codon offset (in nt) of the first
    ATG found, or -1 if none exists in `seq`.
    """
    for i in range(search_init, len(seq) - 2, 3):
        if seq[i:i + 3].upper() == "ATG":
            return i
    return -1


def _compute_next_in_frame(
    annos: pl.LazyFrame,
) -> pl.LazyFrame:
    """Add next_in_frame_relative for start_lost SNVs (polars-native + pyfaidx).

    Only runs on rows where consequence_start_lost == 1. For each such SNV,
    extracts the CDS sequence from the reference FASTA (using cds_start to
    anchor the CDS start on the chromosome, correcting for strand), then
    searches downstream codons for the next ATG. The result is expressed as a
    fraction of CDS length (0–1); variants where no downstream ATG exists get 1.

    Required columns (all from Hail extraction):
        id, chrom, pos, ref, alt, region, cds_start, cds_end, strand
        consequence_start_lost (from step0a)

    cds_start / cds_end semantics (VEP transcript_consequences):
        cds_start  — 1-based nucleotide position of the variant IN the CDS
        cds_end    — same for the last affected nucleotide (== cds_start for SNVs)
        These are NOT the genomic coordinates of the CDS boundaries.

    cds_length is derived from the GTF (sum of CDS exon segments), which is the
    only reliable source. cds_end alone cannot give total CDS length for SNVs
    because it equals cds_start.
    """
    try:
        import pyfaidx  # noqa: PLC0415
    except ImportError:
        logger.warning(
            "  pyfaidx not installed; skipping next_in_frame_relative. "
            "Add pyfaidx to your conda environment."
        )
        return annos

    schema = set(annos.collect_schema().names())
    req = {"id", "chrom", "pos", "ref", "alt", "region",
           "cds_start", "strand", "consequence_start_lost"}
    missing = req - schema
    if missing:
        logger.warning(f"  next_in_frame: missing columns {missing}; skipping")
        return annos

    if "consequence_start_lost" not in schema:
        logger.info("  next_in_frame: consequence_start_lost column absent; skipping")
        return annos

    # ── Subset to start_lost SNVs with usable cds_start ──────────────────
    start_lost = (
        annos
        .filter(pl.col("consequence_start_lost") == 1)
        .filter(
            (pl.col("ref").str.len_chars() == 1) &
            (pl.col("alt").str.len_chars() == 1)
        )
        .filter(pl.col("cds_start").is_not_null())
        .filter(pl.col("strand").is_not_null())
        .select(["id", "region", "chrom", "pos", "cds_start", "strand"])
        .collect()
    )

    if start_lost.is_empty():
        logger.info("  next_in_frame: no start_lost SNVs found; skipping")
        return annos

    # ── Load per-transcript CDS lengths from the GTF ──────────────────────
    # We need CDS length for the denominator. cds_start alone gives the
    # variant's position in the CDS, not the total length.
    gtf_path = os.path.join(WORK_DIR, "gencode.v40.annotation.gtf.gz")
    if not os.path.exists(gtf_path):
        logger.warning("  next_in_frame: GTF not found; skipping")
        return annos

    gtf = _load_gtf_polars(gtf_path)
    # Sum CDS exon lengths per gene (gene-level approximation; transcript-level
    # would need transcript_id join which we don't always have cleanly).
    gene_cds_len = (
        gtf.filter(pl.col("feature") == "CDS")
        .with_columns(
            region=pl.col("attributes")
                   .str.extract(r'gene_id "([^"]+)"')
                   .str.split(".").list.first(),
            seg_len=pl.col("end") - pl.col("start") + 1,
        )
        .group_by("region")
        .agg(cds_length=pl.col("seg_len").sum())
    )
    start_lost = start_lost.join(gene_cds_len, on="region", how="left")

    # ── Fetch sequences and find next ATG ────────────────────────────────
    fasta_gz = os.path.join(WORK_DIR, "GRCh38.primary_assembly.genome.fa.gz")
    if not os.path.exists(fasta_gz):
        if not aria2c_download(GENCODE_FASTA_URL, fasta_gz):
            logger.warning("  next_in_frame: FASTA download failed; skipping")
            return annos
    fasta = pyfaidx.Fasta(fasta_gz, sequence_always_upper=True)
    records = []
    for row in start_lost.iter_rows(named=True):
        var_id    = row["id"]
        chrom     = row["chrom"]
        pos       = row["pos"]           # 1-based genomic position
        cds_pos   = row["cds_start"]     # 1-based position of variant IN CDS
        strand    = str(row["strand"])   # "1" or "-1"
        cds_len   = row.get("cds_length") or 0

        if cds_len == 0:
            records.append({"id": var_id, "region": row["region"],
                             "next_in_frame_relative": 1.0})
            continue

        try:
            chrom_key = chrom if chrom in fasta else chrom.lstrip("chr")
            if strand == "1":
                # Genomic start of CDS = variant_genomic_pos - (cds_pos - 1)
                cds_genome_start = pos - (cds_pos - 1)   # 1-based
                # Fetch CDS + 100 nt buffer for downstream search
                seq = str(fasta[chrom_key][cds_genome_start - 1 : cds_genome_start - 1 + cds_len + 100])
            else:
                # On minus strand VEP cds_pos counts from the transcript 5' end.
                # Genomic end of CDS (highest coordinate) = pos + (cds_pos - 1)
                cds_genome_end = pos + (cds_pos - 1)     # 1-based inclusive
                raw = str(fasta[chrom_key][cds_genome_end - cds_len - 100 : cds_genome_end])
                seq = raw[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))

            dist = _next_inframe_atg_distance(seq, search_init=3)
            rel  = 1.0 if dist < 0 else min(dist / cds_len, 1.0)
        except Exception as fe:
            logger.debug(f"  next_in_frame: FASTA error for {var_id}: {fe}")
            rel = 1.0

        records.append({"id": var_id, "region": row["region"],
                        "next_in_frame_relative": rel})

    if not records:
        return annos

    nif = pl.DataFrame(records, schema={
        "id": pl.Utf8, "region": pl.Utf8, "next_in_frame_relative": pl.Float32
    })
    annos = _idempotent_join(annos, nif.lazy(), on=["id", "region"])
    logger.info(f"  next_in_frame_relative OK ({len(records)} start_lost variants)")
    return annos


def step0b_structural_features(
    annos: pl.LazyFrame,
    work_dir: str,
    add_next_in_frame: bool = True,
) -> pl.LazyFrame:
    """VEP structural features: indel flags, relative_cds_position, dist_to_tss,
    gene_length, gene_name, and (optionally) next_in_frame_relative.

    variant_length / is_indel / is_insertion / is_deletion: from ref/alt lengths.

    relative_cds_position: cds_start / cds_length_from_GTF.
        cds_start (from Hail) = 1-based nt position of the variant within the CDS.
        NOT cds_start / cds_end — for SNVs cds_end == cds_start, making that
        ratio always 1. Total CDS length comes from summing CDS exon segments in
        the Gencode GTF, which is also needed for next_in_frame_relative.

    dist_to_tss / gene_length / gene_name: from the Gencode GTF, joined on region.

    next_in_frame_relative: downloaded and computed automatically for start_lost SNVs.
        Pass add_next_in_frame=False to main() to skip entirely.
    """
    logger.info("Step 0b: VEP structural features")
    schema = set(annos.collect_schema().names())

    # ── variant length + indel flags ──────────────────────────────────────
    annos = annos.with_columns(
        variant_length=pl.max_horizontal(
            pl.col("ref").str.len_chars(), pl.col("alt").str.len_chars()
        )
    ).with_columns(
        is_indel=(pl.col("variant_length") > 1).cast(pl.Int8),
        is_insertion=(
            pl.col("ref").str.len_chars() < pl.col("alt").str.len_chars()
        ).cast(pl.Int8),
        is_deletion=(
            pl.col("ref").str.len_chars() > pl.col("alt").str.len_chars()
        ).cast(pl.Int8),
    )

    # ── GTF: needed for both relative_cds_position and next_in_frame ──────
    gtf_path = os.path.join(work_dir, "gencode.v40.annotation.gtf.gz")
    if not os.path.exists(gtf_path):
        if not aria2c_download(GENCODE_GTF_URL, gtf_path):
            logger.warning("  Gencode GTF download failed; skipping GTF-derived features")
            return annos

    try:
        gtf = _load_gtf_polars(gtf_path)

        # ── dist_to_tss, gene_length, gene_name ───────────────────────────
        genes = _gtf_gene_features(gtf)
        annos = (
            _idempotent_join(annos, genes.lazy(), on=["region"], validate="m:1")
            .with_columns(
                dist_to_tss=pl.when(pl.col("_gene_strand") == "+")
                .then(pl.col("pos") - pl.col("_tss"))
                .otherwise(pl.col("_tss") - pl.col("pos"))
            )
            .drop(["_tss", "_gene_strand"])
        )
        logger.info("  dist_to_tss / gene_length / gene_name OK")

        # ── relative_cds_position: cds_start / total CDS length ───────────
        # cds_start is the 1-based nt position of the variant within the CDS
        # (from VEP transcript_consequences). Total CDS length = sum of CDS
        # exon segments from the GTF (gene-level; avoids needing transcript_id).
        if "cds_start" in schema:
            gene_cds_len = (
                gtf.filter(pl.col("feature") == "CDS")
                .with_columns(
                    region=pl.col("attributes")
                           .str.extract(r'gene_id "([^"]+)"')
                           .str.split(".").list.first(),
                    seg_len=pl.col("end") - pl.col("start") + 1,
                )
                .group_by("region")
                .agg(cds_length=pl.col("seg_len").sum())
            )
            annos = (
                _idempotent_join(annos, gene_cds_len.lazy(), on=["region"])
                .with_columns(
                    relative_cds_position=(
                        pl.col("cds_start").cast(pl.Float64)
                        / pl.col("cds_length").cast(pl.Float64)
                    ).clip(0.0, 1.0).round(4).cast(pl.Float32)
                )
                # keep cds_length — also used by next_in_frame
            )
            logger.info("  relative_cds_position OK (cds_start / GTF CDS length)")
        else:
            logger.warning(
                "  cds_start absent; skipping relative_cds_position. "
                "Extract from vep.transcript_consequences[].cds_start in Hail."
            )
    except Exception as e:
        logger.warning(f"  GTF-derived features failed: {e}")

    # ── next_in_frame_relative (start_lost SNVs only) ────────────────────
    if add_next_in_frame:
        annos = _compute_next_in_frame(annos)

    return annos


# ── Main entry point ───────────────────────────────────────────────────────
def main(
    vep_parquet: str,
    fill_null_defaults_path: str = "fill_null_defaults.yaml",
    output_path: str | None = None,
    download_dir: str | None = None,
    *,
    add_alphamissense: bool = True,
    add_popeve: bool = True,
    add_revel: bool = True,
    add_clinpred: bool = True,
    add_bayesdel: bool = True,
    add_cpt1: bool = True,
    add_cadd: bool = True,
    add_gpn_msa: bool = True,
    add_phylop: bool = True,
    add_plddt: bool = True,
    add_pioneer: bool = True,
    add_clinvar: bool = True,
    add_next_in_frame: bool = True,
) -> str:
    """Annotate the VEP parquet with effect-prediction scores and write it out.

    Parameters
    ----------
    vep_parquet:
        Path (local or gs://) to the vep_loftee_parallel output parquet.
        Read directly with polars (gs:// requires gcsfs to be installed).
    fill_null_defaults_path:
        Path to fill_null_defaults.yaml (shipped alongside this script in the
        repo). Passed explicitly because __file__ is undefined in notebooks.
    output_path:
        Where to write the final parquet. Defaults to
        ./{input_stem}_annotated.parquet in the current working directory.
    download_dir:
        Directory under which the large reference/score files are downloaded.
        They land in ``{download_dir}/tmp``. That ``tmp`` dir is deleted only
        after a fully successful run; if the pipeline errors out it is kept so
        the rerun reuses already-downloaded files. Defaults to the current
        working directory, i.e. ``./tmp``.
    add_alphamissense, add_revel, add_clinpred, add_bayesdel, add_cpt1,
    add_cadd, add_gpn_msa, add_phylop, add_plddt, add_pioneer, add_clinvar:
        Per-annotation toggles (keyword-only, all default True). Set one to
        False to skip both that annotation's download and its merge step.
        Useful when the input parquet already carries a column (e.g. the
        DeepRVAT WGS snakemake output already has am_pathogenicity, cpt1_llr,
        cadd_*, gpn_score and the phyloP/conservation columns): re-running those
        steps is redundant (an idempotent join would keep the existing column
        anyway) and only wastes the multi-GB CADD/GPN-MSA/bigWig downloads and
        the slow tabix/bigWig queries. The genuinely new annotations this script
        adds on top of that schema are REVEL, ClinPred, BayesDel, pLDDT, PIONEER
        and ClinVar. ``add_phylop`` controls all three conservation bigWigs
        (vertebrate 100-way, mammalian 17-way, primate 30-way) together.

    Returns
    -------
    The path to the written parquet.
    """
    # Route all downloads into a scratch tmp dir. WORK_DIR is a module global
    # referenced by the download helpers and step functions, so rebind it here
    # for the duration of this run. It is deleted only on success (see end) —
    # on failure the cache survives so reruns skip what's already downloaded.
    global WORK_DIR
    base_dir = download_dir if download_dir is not None else os.getcwd()
    WORK_DIR = os.path.join(base_dir, "tmp")
    os.makedirs(WORK_DIR, exist_ok=True)

    # ── Load fill-null defaults ──────────────────────────────────────
    with open(fill_null_defaults_path) as f:
        fill_defaults = yaml.safe_load(f)

    # ── Required downloads (always attempted) ───────────────────────────
    logger.info("=== Downloading required files ===")

    # On AoU the VEP parquet is read directly from its path (local or gs://);
    # there is no dxpy download step. polars reads gs:// when gcsfs is installed.
    vep_path = vep_parquet

    am_path = f"{WORK_DIR}/AlphaMissense_hg38.tsv.gz"
    am_ok = add_alphamissense and aria2c_download(AM_URL, am_path)

    clinvar_path = f"{WORK_DIR}/clinvar.vcf.gz"
    clinvar_ok = add_clinvar and aria2c_download(CLINVAR_URL, clinvar_path)

    # ── Auto-download annotation scores ─────────────────────────────────
    logger.info("=== Auto-downloading annotation scores ===")

    # REVEL (hg38 CSV, ~1.4 GB zip)
    revel_path = None
    if add_revel:
        revel_zip = f"{WORK_DIR}/revel.zip"
        if aria2c_download(REVEL_URL, revel_zip):
            revel_dir = unzip(revel_zip, f"{WORK_DIR}/revel")
            # The Zenodo zip contains a single file named "revel_with_transcript_ids"
            # (no extension). Fall back to *.csv for older mirrors.
            candidates = (
                list(Path(revel_dir).rglob("revel_with_transcript_ids"))
                + list(Path(revel_dir).rglob("*.csv"))
            )
            revel_path = str(candidates[0]) if candidates else None

    # ClinPred (hg38 coordinates)
    clinpred_path = None
    if add_clinpred:
        clinpred_gz = f"{WORK_DIR}/ClinPred_hg38.txt.gz"
        if gdrive_download(CLINPRED_GDRIVE_ID, clinpred_gz):
            clinpred_path = clinpred_gz

    # BayesDel noAF — downloaded as .gz; polars scan_csv reads gzip natively
    bayesdel_path = None
    if add_bayesdel:
        bayesdel_gz = f"{WORK_DIR}/bayesdel.gz"
        if gdrive_download(BAYESDEL_GDRIVE_ID, bayesdel_gz):
            bayesdel_path = bayesdel_gz

    # popEVE / EVE / ESM1v per-variant VCF (~1.4 GB)
    popeve_path = None
    if add_popeve:
        popeve_gz = f"{WORK_DIR}/grch38_popEVE.vcf.gz"
        if aria2c_download(POPEVE_URL, popeve_gz):
            popeve_path = popeve_gz

    # CPT-1 per-protein CSVs from Zenodo 8140323 (three zips, ~2.4 GB total)
    cpt1_dir = f"{WORK_DIR}/cpt1"
    cpt1_ok = False
    if add_cpt1:
        os.makedirs(cpt1_dir, exist_ok=True)
        for idx, url in enumerate(CPT1_ZENODO_URLS):
            zip_path = f"{WORK_DIR}/cpt1_{idx}.zip"
            if aria2c_download(url, zip_path):
                unzip(zip_path, cpt1_dir)
                cpt1_ok = True

    # UniProt ID mapping (used by CPT-1 and the protein-structure steps).
    # Only needed when at least one of those steps will run.
    if add_cpt1 or add_plddt or add_pioneer:
        logger.info("=== Building UniProt map ===")
        uniprot_map_df = _build_uniprot_map(WORK_DIR)
    else:
        uniprot_map_df = pl.DataFrame(
            schema={"region": pl.Utf8, "uniprot_id": pl.Utf8,
                    "gene_name": pl.Utf8, "entry_name": pl.Utf8}
        )

    # NOTE: Splicing predictors (Pangolin, AbSplice, AbSplice2), AbExp expression
    # scores and PromoterAI are out of scope for this exome/missense effector
    # port and have been removed. Their download blocks and step functions are
    # gone; re-add them when you extend to the non-coding / splicing tracks.

    # ── Load VEP parquet, rename gene_id → region ────────────────────────
    logger.info("=== Loading VEP parquet ===")
    annos = pl.scan_parquet(vep_path)
    schema_names = set(annos.collect_schema().names())
    if "gene_id" in schema_names and "region" not in schema_names:
        annos = annos.rename({"gene_id": "region"})

    # ── Harmonize AoU schema to what the downstream steps expect ─────────
    # - pos → Int64 (AoU ships Int32; the score-table joins use Int64 keys)
    # - chrom → 'chr'-prefixed (idempotent), matching every join below
    # - id → 'chrom:pos:ref:alt' (same format step4d_clinvar builds; required
    #   because the AoU file has variant_id, not id, and several steps select 'id')
    # - protein_position → str(protein_start) so CPT-1 / pLDDT / PIONEER, which
    #   parse 'protein_position', resolve a residue number on AoU rows
    logger.info("=== Harmonizing AoU schema ===")
    annos = annos.with_columns(
        pl.col("pos").cast(pl.Int64),
        chrom=pl.lit("chr") + pl.col("chrom").cast(pl.Utf8).str.replace(r"^chr", ""),
    )
    schema_names = set(annos.collect_schema().names())
    if "id" not in schema_names:
        annos = annos.with_columns(
            id=pl.col("chrom") + ":" + pl.col("pos").cast(pl.Utf8)
            + ":" + pl.col("ref") + ":" + pl.col("alt")
        )
    if "protein_position" not in schema_names and "protein_start" in schema_names:
        annos = annos.with_columns(
            protein_position=pl.col("protein_start").cast(pl.Utf8)
        )

    # ── Attach absolute protein_position / amino_acids from raw VEP ──────
    n_rows = annos.select(pl.len()).collect().item()
    logger.info(f"Loaded {n_rows:,} rows")

    # ── Group 0: VEP / LOFTEE post-processing ────────────────────────────
    logger.info("=== Group 0: VEP / LOFTEE post-processing ===")
    annos = step0a_vep_loftee_dummies(annos)
    annos = step0b_structural_features(annos, WORK_DIR, add_next_in_frame=add_next_in_frame)

    # ── Group 1: Missense / protein-function predictors ─────────────────
    logger.info("=== Group 1: Missense predictors ===")
    if am_ok:
        annos = step1a_alphamissense(annos, am_path, WORK_DIR)
    if revel_path:
        annos = step1b_revel(annos, revel_path, WORK_DIR)
    if clinpred_path:
        annos = step1c_clinpred(annos, clinpred_path)
    if bayesdel_path:
        annos = step1d_bayesdel(annos, bayesdel_path)
    if add_popeve and popeve_path:
        annos = step1e_popeve(annos, popeve_path)
    if cpt1_ok:
        annos = step1f_cpt1(annos, cpt1_dir, uniprot_map_df)

    # ── Group 2: CADD / conservation ────────────────────────────────────
    logger.info("=== Group 2: CADD / conservation ===")
    if add_cadd:
        annos = step2a_cadd(annos)
    if add_gpn_msa:
        annos = step2b_gpn_msa(annos)
    if add_phylop:
        annos = step2c1_phylop(annos)
        annos = step2c2_phylop_mammalian(annos)
        annos = step2c3_phylop_primate(annos)

    # ── Group 3: Derived columns + fill nulls ───────────────────────────
    logger.info("=== Group 3: Derived columns + fill nulls ===")
    # ENCODE cCRE annotation (step4c_encode) is regulatory/non-coding and is out
    # of scope for the exome port, so it is NOT called here. The function is left
    # defined for later use. To re-enable: define CCRE_COLUMN_MAP, CCRE_ZERO_COLS_V4
    # and ENCODE_URL, download the cCRE BED, and call step4c_encode BEFORE
    # step3a_derived_columns (which reads the encode_* columns it produces).
    annos = step3a_derived_columns(annos)
    annos = step3b_fill_nulls(annos, fill_defaults)

    # ── Group 4: Post-fill annotations (no null handling needed) ────────
    logger.info("=== Group 4: Post-fill annotations ===")
    if add_plddt:
        annos = step4a_plddt(annos, uniprot_map_df, WORK_DIR)
    if add_pioneer:
        annos = step4b_pioneer_interface(annos, uniprot_map_df, WORK_DIR)
    if clinvar_ok:
        annos = step4d_clinvar(annos, clinvar_path)

    # ── Collect and write ─────────────────────────────────────────────────
    logger.info("=== Sinking output parquet ===")
    gc.collect()
    if output_path is None:
        stem = os.path.splitext(os.path.basename(vep_parquet))[0]
        output_path = f"{stem}_annotated.parquet"
    annos.sink_parquet(output_path)
    logger.info(f"Output written to {output_path}")

    # Success: remove the scratch download dir. (On failure we never reach
    # here, so the cache survives for a fast rerun.)
    shutil.rmtree(WORK_DIR, ignore_errors=True)
    logger.info(f"Removed scratch download dir {WORK_DIR}")
    return output_path


# In a notebook, call directly, e.g.:
#   out = main("gs://.../vep_annotated.parquet",
#              fill_null_defaults_path="fill_null_defaults.yaml")