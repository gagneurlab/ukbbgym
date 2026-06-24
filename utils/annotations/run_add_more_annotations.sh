#!/bin/bash
#SBATCH --job-name=ukbgym_annos
#SBATCH --output=/s/project/ukbbgym/annotation_files/more_annotations/logs/%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00

# ── Usage ─────────────────────────────────────────────────────────────────────
#
#   sbatch run_add_more_annotations.sh --input  /s/project/deeprvat_wgs/input_data/annotations/qced_maf1e-3_loftee_olink_genes_EURunrelated/annotations_with_all_no_dup.parquet --output /s/project/deeprvat_wgs/input_data/annotations/qced_maf1e-3_loftee_olink_genes_EURunrelated/annotations_no_dup_20260624.parquet --download-dir /s/project/ukbbgym/annotation_files/more_annotations --no-cadd --no-gpn-msa --no-phylop --no-alphamissense --no-cpt1 --no-next-in-frame
#
#   Or run directly — activate the correct conda env first, then:
#       bash run_add_more_annotations.sh --input ... --output ...
#
# Per-annotation toggles (pass to skip steps whose columns are already present):
#   --no-cadd --no-gpn-msa --no-phylop --no-alphamissense --no-cpt1
#   --no-popeve --no-revel --no-clinpred --no-bayesdel
#   --no-plddt --no-pioneer --no-clinvar --no-next-in-frame
#
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

if [[ -n "${SLURM_JOB_ID:-}" && -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

# ── Defaults ──────────────────────────────────────────────────────────────────
INPUT=""
OUTPUT=""
DOWNLOAD_DIR="/s/project/ukbbgym/annotation_files/more_annotations"

ADD_CADD=True
ADD_GPN_MSA=True
ADD_PHYLOP=True
ADD_ALPHAMISSENSE=True
ADD_CPT1=True
ADD_POPEVE=True
ADD_REVEL=True
ADD_CLINPRED=True
ADD_BAYESDEL=True
ADD_PLDDT=True
ADD_PIONEER=True
ADD_CLINVAR=True
ADD_NEXT_IN_FRAME=True

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)          INPUT="$2";        shift 2 ;;
        --output)         OUTPUT="$2";       shift 2 ;;
        --download-dir)   DOWNLOAD_DIR="$2"; shift 2 ;;
        --no-alphamissense) ADD_ALPHAMISSENSE=False; shift ;;
        --no-popeve)        ADD_POPEVE=False;        shift ;;
        --no-revel)         ADD_REVEL=False;         shift ;;
        --no-clinpred)      ADD_CLINPRED=False;      shift ;;
        --no-bayesdel)      ADD_BAYESDEL=False;      shift ;;
        --no-cpt1)          ADD_CPT1=False;          shift ;;
        --no-cadd)          ADD_CADD=False;          shift ;;
        --no-gpn-msa)       ADD_GPN_MSA=False;       shift ;;
        --no-phylop)        ADD_PHYLOP=False;         shift ;;
        --no-plddt)         ADD_PLDDT=False;         shift ;;
        --no-pioneer)       ADD_PIONEER=False;       shift ;;
        --no-clinvar)       ADD_CLINVAR=False;        shift ;;
        --no-next-in-frame) ADD_NEXT_IN_FRAME=False;  shift ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$INPUT" ]]; then
    echo "Error: --input is required" >&2
    exit 1
fi

# ── Build optional keyword args ───────────────────────────────────────────────
OPTIONAL_ARGS=""
[[ -n "$OUTPUT"       ]] && OPTIONAL_ARGS+="    output_path='${OUTPUT}',"$'\n'
[[ -n "$DOWNLOAD_DIR" ]] && OPTIONAL_ARGS+="    download_dir='${DOWNLOAD_DIR}',"$'\n'

echo "=== add_more_annotations ==="
echo "  input:        $INPUT"
echo "  output:       ${OUTPUT:-<auto>}"
echo "  download_dir: ${DOWNLOAD_DIR}"
echo "  toggles: alphamissense=$ADD_ALPHAMISSENSE popeve=$ADD_POPEVE revel=$ADD_REVEL"
echo "           clinpred=$ADD_CLINPRED bayesdel=$ADD_BAYESDEL cpt1=$ADD_CPT1"
echo "           cadd=$ADD_CADD gpn_msa=$ADD_GPN_MSA phylop=$ADD_PHYLOP"
echo "           plddt=$ADD_PLDDT pioneer=$ADD_PIONEER clinvar=$ADD_CLINVAR next_in_frame=$ADD_NEXT_IN_FRAME"
echo ""

mkdir -p /s/project/ukbbgym/annotation_files/more_annotations/logs

python3 - <<PYEOF
import sys
sys.path.insert(0, '${SCRIPT_DIR}')
import add_more_annotations as ann

out = ann.main(
    '${INPUT}',
    fill_null_defaults_path='${SCRIPT_DIR}/fill_null_defaults.yaml',
${OPTIONAL_ARGS}    add_alphamissense=${ADD_ALPHAMISSENSE},
    add_popeve=${ADD_POPEVE},
    add_revel=${ADD_REVEL},
    add_clinpred=${ADD_CLINPRED},
    add_bayesdel=${ADD_BAYESDEL},
    add_cpt1=${ADD_CPT1},
    add_cadd=${ADD_CADD},
    add_gpn_msa=${ADD_GPN_MSA},
    add_phylop=${ADD_PHYLOP},
    add_plddt=${ADD_PLDDT},
    add_pioneer=${ADD_PIONEER},
    add_clinvar=${ADD_CLINVAR},
    add_next_in_frame=${ADD_NEXT_IN_FRAME},
)
print(f"Done: {out}")
PYEOF
