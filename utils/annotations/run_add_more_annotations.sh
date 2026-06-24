#!/bin/bash
#SBATCH --job-name=add_missense_annotations
#SBATCH --output=logs/add_missense_annotations/%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=cpu_p
#SBATCH --time=24:00:00

# ── Usage ─────────────────────────────────────────────────────────────────────
#
#   sbatch run_add_missense_annotations.sh \
#       --input  /path/to/vep.parquet \
#       --output /path/to/output_annotated.parquet \
#       [--download-dir /path/to/download_cache] \
#       [--vep-raw   /path/to/raw_vep.parquet]
#
#   Or run directly (no SLURM):
#       bash run_add_missense_annotations.sh --input ... --output ...
#
# Per-annotation toggles (pass to skip a step whose columns are already present):
#   --no-cadd --no-gpn-msa --no-phylop --no-alphamissense --no-cpt1
#   --no-popeve --no-revel --no-clinpred --no-bayesdel
#   --no-plddt --no-pioneer --no-clinvar
#
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON=/opt/modules/i12g/anaconda/envs/sl-ukg/bin/python3

# ── Defaults ──────────────────────────────────────────────────────────────────
INPUT=""
OUTPUT=""
DOWNLOAD_DIR=""
VEP_RAW=""
FASTA=""

ADD_ALPHAMISSENSE=False
ADD_CPT1=False
ADD_CADD=False
ADD_GPN_MSA=False
ADD_PHYLOP=False
ADD_POPEVE=True
ADD_REVEL=True
ADD_CLINPRED=True
ADD_BAYESDEL=True
ADD_PLDDT=True
ADD_PIONEER=True
ADD_CLINVAR=True

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)          INPUT="$2";        shift 2 ;;
        --output)         OUTPUT="$2";       shift 2 ;;
        --download-dir)   DOWNLOAD_DIR="$2"; shift 2 ;;
        --vep-raw)        VEP_RAW="$2";      shift 2 ;;
        --fasta)          FASTA="$2";        shift 2 ;;
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
[[ -n "$VEP_RAW"      ]] && OPTIONAL_ARGS+="    vep_raw_parquet='${VEP_RAW}',"$'\n'
[[ -n "$FASTA"        ]] && OPTIONAL_ARGS+="    fasta_path='${FASTA}',"$'\n'

echo "=== add_missense_variant_annotations ==="
echo "  input:        $INPUT"
echo "  output:       ${OUTPUT:-<auto>}"
echo "  download_dir: ${DOWNLOAD_DIR:-<cwd>/tmp}"
echo "  vep_raw:      ${VEP_RAW:-<none>}"
echo "  toggles: alphamissense=$ADD_ALPHAMISSENSE popeve=$ADD_POPEVE revel=$ADD_REVEL"
echo "           clinpred=$ADD_CLINPRED bayesdel=$ADD_BAYESDEL cpt1=$ADD_CPT1"
echo "           cadd=$ADD_CADD gpn_msa=$ADD_GPN_MSA phylop=$ADD_PHYLOP"
echo "           plddt=$ADD_PLDDT pioneer=$ADD_PIONEER clinvar=$ADD_CLINVAR"
echo ""

mkdir -p logs/add_missense_annotations

"$PYTHON" - <<PYEOF
import sys
sys.path.insert(0, '${SCRIPT_DIR}')
import add_missense_variant_annotations as ann

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
)
print(f"Done: {out}")
PYEOF
