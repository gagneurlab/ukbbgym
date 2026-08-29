#!/usr/bin/env python3
"""CLI used by run_all.sh --new-score: merge a new predictor column onto the master table
(merge_new_score) and register it in config_correlations.yaml (register_score_in_config),
both in utils/variant_filtering.py. Prints MASTER_PATH=... and CONFIG_FILE=... for the
caller to export.
"""

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
from utils.variant_filtering import fetch_hf_data, merge_new_score, register_score_in_config


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('score_path')
    ap.add_argument('--master-path', default=None)
    ap.add_argument('--category', default='missense')
    ap.add_argument('--on', default=None, help='force the join key (name, or comma-separated '
                     'for a compound key) instead of auto-detecting genomic vs. protein')
    args = ap.parse_args()

    master_path = args.master_path or fetch_hf_data('genebass_annotated.parquet', REPO_ROOT)
    on = args.on.split(',') if args.on else None
    merged_path, new_cols = merge_new_score(master_path, args.score_path, on=on)

    cfg_name = register_score_in_config(str(REPO_ROOT / 'configs'), 'config_correlations.yaml',
                                         new_cols, category=args.category)
    print(f'MASTER_PATH={merged_path}')
    print(f'CONFIG_FILE={cfg_name}')
    print(f"Registered {new_cols} under category '{args.category}' in {cfg_name}", file=sys.stderr)


if __name__ == '__main__':
    main()
