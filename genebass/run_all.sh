#!/usr/bin/env bash
# Run every genebass/analysis notebook top to bottom, in place, regenerating all figures
# in paper_figures/ at the repo root. Requires the master table at data/ (see README.md);
# some notebooks also need the additional data/ inputs listed in the README's inventory --
# those are allowed to fail here (reported at the end) rather than aborting the whole run.
#
# Usage: genebass/run_all.sh
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PYTHON="$REPO_ROOT/.venv/bin/python"

if [[ ! -x "$PYTHON" ]]; then
    echo "No .venv at $REPO_ROOT/.venv -- run 'uv sync' from the repo root first." >&2
    exit 1
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
