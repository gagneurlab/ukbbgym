#!/usr/bin/env bash
# Run every genebass/analysis notebook top to bottom, in place, regenerating all figures
# in paper_figures/ at the repo root. Each notebook auto-fetches its own data/ inputs from
# the private gagneurlab/ukbbgym Hugging Face dataset on first run (see README.md) -- this
# requires a logged-in HF account with access; a 401 there is allowed to fail a notebook
# here (reported at the end) rather than aborting the whole run.
#
# Each notebook's parameter cell reads its defaults via env_override() (utils/variant_filtering.py),
# so any of the flags below overrides that one setting in every notebook that has it, without
# editing the notebooks -- a notebook that has no such parameter just ignores the flag. A few
# notebooks give a differently-named parameter its own UKBBGYM_<NAME> instead of sharing one of
# these flags, specifically because its default means something different there (a different
# config file's categories, or a different master table) -- see genebass/README.md's table for
# those; set them by exporting the env var yourself, there's no flag for them.
#
# Registers a `ukbbgym-venv` Jupyter kernel pointing at this repo's .venv on every run (see
# below) -- the standard fix for nbconvert otherwise picking whatever kernel happens to be
# registered globally; this is the only effect outside the repo.
#
# Usage: genebass/run_all.sh [options]
#   --variant-class NAME        e.g. missense, indel, missense_structured (default: per-notebook)
#   --selected-categories LIST  comma-separated config_correlations.yaml categories, e.g.
#                                "missense,conservation" (default: per-notebook)
#   --mac N                     rare-variant allele-count cap (default: 20)
#   --min-variants N            coverage threshold below which a gene/pair is dropped (default: per-notebook)
#   --config-file NAME          e.g. config_correlations.yaml (default: per-notebook)
#   --master-path PATH          override the shared master table path (default: auto-fetched from Hugging Face into data/genebass_annotated.parquet -- requires a logged-in HF account, see README.md)
#   --fig-dir PATH              where figures are written (default: paper_figures/)
#   --only-snps / --no-only-snps
#   --new-score PATH            merge a new predictor column onto the master table before running
#                                (parquet/csv with an 'id' column plus one or more score columns --
#                                see merge_new_score.py) and register it as a predictor, so every
#                                notebook's plots include it alongside the existing tools. Default:
#                                none -- run with the existing annotations only.
#   --new-score-category NAME   config_correlations.yaml category the new score is registered
#                                under (default: missense -- included in every notebook's default
#                                selected_categories, so it shows up without --selected-categories)
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PYTHON="$REPO_ROOT/.venv/bin/python"

NEW_SCORE=""
NEW_SCORE_CATEGORY="missense"

while (($#)); do
    case "$1" in
        --variant-class)        export UKBBGYM_VARIANT_CLASS="$2"; shift 2 ;;
        --selected-categories)  export UKBBGYM_SELECTED_CATEGORIES="$2"; shift 2 ;;
        --mac)                  export UKBBGYM_MAC="$2"; shift 2 ;;
        --min-variants)         export UKBBGYM_MIN_VARIANTS="$2"; shift 2 ;;
        --config-file)          export UKBBGYM_CONFIG_FILE="$2"; shift 2 ;;
        --master-path)          export UKBBGYM_MASTER_PATH="$2"; shift 2 ;;
        --fig-dir)               export UKBBGYM_FIG_DIR="$2"; shift 2 ;;
        --only-snps)             export UKBBGYM_ONLY_SNPS=true; shift ;;
        --no-only-snps)          export UKBBGYM_ONLY_SNPS=false; shift ;;
        --new-score)             NEW_SCORE="$2"; shift 2 ;;
        --new-score-category)    NEW_SCORE_CATEGORY="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

if [[ ! -x "$PYTHON" ]]; then
    echo "No .venv at $REPO_ROOT/.venv -- run 'uv sync' from the repo root first." >&2
    exit 1
fi

if [[ -n "$NEW_SCORE" ]]; then
    echo "=== merging new score: $NEW_SCORE ==="
    merge_out="$("$PYTHON" "$SCRIPT_DIR/merge_new_score.py" "$NEW_SCORE" \
        --category "$NEW_SCORE_CATEGORY" \
        ${UKBBGYM_MASTER_PATH:+--master-path "$UKBBGYM_MASTER_PATH"})" || exit 1
    export UKBBGYM_MASTER_PATH="$(sed -n 's/^MASTER_PATH=//p' <<< "$merge_out")"
    export UKBBGYM_CONFIG_FILE="$(sed -n 's/^CONFIG_FILE=//p' <<< "$merge_out")"
    echo "  master table -> $UKBBGYM_MASTER_PATH"
    echo "  config file  -> $UKBBGYM_CONFIG_FILE"
fi

# nbconvert picks a kernel by name from Jupyter's global kernelspec list, not from
# whichever python invokes it -- without this it can silently run notebooks against some
# other kernel on the machine. Register (or refresh) a kernel pointing at this repo's .venv.
KERNEL=ukbbgym-venv
"$PYTHON" -m ipykernel install --user --name "$KERNEL" --display-name "ukbbgym (.venv)" >/dev/null

mapfile -t NOTEBOOKS < <(find "$SCRIPT_DIR/analysis" -name '*.ipynb' | sort)

failed=()
for nb in "${NOTEBOOKS[@]}"; do
    rel="${nb#"$REPO_ROOT"/}"
    echo "=== $rel ==="
    if "$PYTHON" -m jupyter nbconvert --to notebook --execute --inplace \
            --ExecutePreprocessor.kernel_name="$KERNEL" \
            --ExecutePreprocessor.timeout=1800 "$nb"; then
        echo "  OK"
    else
        echo "  FAILED"
        failed+=("$rel")
    fi
done

echo
echo "Figures written to $REPO_ROOT/paper_figures/"
if ((${#failed[@]})); then
    echo "${#failed[@]} notebook(s) failed (likely missing data/ inputs -- see genebass/README.md):"
    printf '  %s\n' "${failed[@]}"
    exit 1
fi
echo "All notebooks ran clean."
