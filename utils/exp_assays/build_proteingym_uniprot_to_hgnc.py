"""Rebuild proteingym_uniprot_to_hgnc.csv from the UniProt REST API.

ProteinGym DMS substitution files are named like ``HXK4_HUMAN_Gersing_2022_activity.csv``.
The pipeline keys on the first underscore-delimited token (``HXK4``), which is the prefix of
a UniProt *entry name* (``HXK4_HUMAN``) -- not an HGNC gene symbol (``GCK``). About 45% of the
human assays disagree. This script resolves every token to its primary gene symbol by asking
UniProt, so the mapping is reproducible rather than hand-maintained.

Usage:
    python build_proteingym_uniprot_to_hgnc.py \
        --dms-dir /path/to/DMS_ProteinGym_substitutions \
        --out proteingym_uniprot_to_hgnc.csv

Needs network access to https://rest.uniprot.org.
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"

# Tokens UniProt cannot resolve from "<token>_HUMAN" alone (e.g. it is an accession, or the
# entry name changed). Checked into the script so a rerun is deterministic; extend as needed.
FALLBACK = {
    "Q53Z42": "HLA-A",  # token is a UniProt accession, not an entry-name prefix
}


def tokens_from_dms_dir(dms_dir: Path) -> list[str]:
    toks = {
        f.name.split(".")[0].split("_")[0]
        for f in dms_dir.glob("*.csv")
        if "HUMAN" in f.name
    }
    return sorted(toks)


def query_uniprot(token: str) -> str | None:
    """Primary gene symbol for entry name '<token>_HUMAN' (or accession '<token>'), or None."""
    query = f"(id:{token}_HUMAN) OR (accession:{token})"
    url = (
        f"{UNIPROT_SEARCH}?"
        + urllib.parse.urlencode(
            {"query": query, "fields": "gene_primary", "format": "tsv", "size": "1"}
        )
    )
    with urllib.request.urlopen(url, timeout=30) as resp:
        rows = list(csv.reader(io.StringIO(resp.read().decode()), delimiter="\t"))
    if len(rows) < 2 or not rows[1] or not rows[1][0].strip():
        return None
    # gene_primary can be "GENE; SYNONYM"; take the first.
    return rows[1][0].split(";")[0].strip()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dms-dir", type=Path, required=True,
                    help="ProteinGym DMS_ProteinGym_substitutions directory")
    ap.add_argument("--out", type=Path, default=Path("proteingym_uniprot_to_hgnc.csv"))
    ap.add_argument("--sleep", type=float, default=0.2, help="seconds between API calls")
    args = ap.parse_args()

    toks = tokens_from_dms_dir(args.dms_dir)
    if not toks:
        sys.exit(f"no *HUMAN*.csv files under {args.dms_dir}")
    print(f"{len(toks)} filename tokens", file=sys.stderr)

    rows, unresolved = [], []
    for t in toks:
        sym = FALLBACK.get(t) or query_uniprot(t)
        if sym is None:
            unresolved.append(t)
            sym = t  # leave the token as-is; flagged below
        rows.append((t, sym))
        print(f"  {t:8s} -> {sym}", file=sys.stderr)
        time.sleep(args.sleep)

    with args.out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["filename_token", "hgnc_symbol"])
        w.writerows(rows)

    print(f"wrote {len(rows)} rows -> {args.out}", file=sys.stderr)
    if unresolved:
        print(f"UNRESOLVED (left as token, add to FALLBACK): {unresolved}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
