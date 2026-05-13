#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="$(dirname "$(dirname "${BASH_SOURCE[0]}")")/share/molsimstack"
export PYTHONPATH="$BASE_DIR:${PYTHONPATH:-}"

INPUT=""
OUTDIR=""
METHODS=()
THREADS=1
DOMAIN=""

usage() {
    echo "MolSimStack - Molecular Similarity Pipeline"
    echo ""
    echo "Usage:"
    echo "  molsimstack.sh -i <input> -o <outdir> -d <domain> [-m methods] [-t threads]"
    echo ""
    echo "Domains:"
    echo "  nuc   → nucleic acids (BLAST-based)"
    echo "  mol   → small molecules / peptides (RDKit-based)"
    echo ""
    exit 1
}

[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            [[ $# -lt 2 ]] && { echo "ERROR: Missing value for $1"; exit 1; }
            INPUT="$2"; shift 2 ;;
        -o|--outdir)
            [[ $# -lt 2 ]] && { echo "ERROR: Missing value for $1"; exit 1; }
            OUTDIR="$2"; shift 2 ;;
        -d|--domain)
            [[ $# -lt 2 ]] && { echo "ERROR: Missing value for $1"; exit 1; }
            DOMAIN="$2"; shift 2 ;;
        -m|--methods)
            shift
            while [[ $# -gt 0 && ! "$1" =~ ^- ]]; do
                METHODS+=("$1")
                shift
            done ;;
        -t|--threads)
            [[ $# -lt 2 ]] && { echo "ERROR: Missing value for $1"; exit 1; }
            THREADS="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "ERROR: Unknown argument: $1"; usage ;;
    esac
done

# -----------------------------
# Required arguments
# -----------------------------
[[ -z "$INPUT" || -z "$OUTDIR" || -z "$DOMAIN" ]] && usage

# -----------------------------
# Input path validation (file OR directory)
# -----------------------------
if [[ ! -f "$INPUT" && ! -d "$INPUT" ]]; then
    echo "ERROR: Input path not found (must be file or directory)"
    exit 1
fi

# -----------------------------
# Domain validation
# -----------------------------
if [[ "$DOMAIN" != "nuc" && "$DOMAIN" != "mol" ]]; then
    echo "ERROR: domain must be 'nuc' or 'mol'"
    exit 1
fi

# -----------------------------
# Molecular input validation
# -----------------------------
if [[ "$DOMAIN" == "mol" ]]; then
    # If it's a file → must be .sdf
    if [[ -f "$INPUT" && "$INPUT" != *.sdf ]]; then
        echo "ERROR: molecular domain requires .sdf file or directory of .sdf files"
        exit 1
    fi
fi

# -----------------------------
# Nucleic input validation
# -----------------------------
if [[ "$DOMAIN" == "nuc" ]]; then
    # If it's a file → must be FASTA
    if [[ -f "$INPUT" && "$INPUT" != *.fasta && "$INPUT" != *.fa ]]; then
        echo "ERROR: nucleic domain requires FASTA file or directory"
        exit 1
    fi
fi


# BLAST check
if [[ "$DOMAIN" == "nuc" ]]; then
    command -v blastn >/dev/null || { echo "ERROR: blastn not found"; exit 1; }
fi

mkdir -p "$OUTDIR"

METHOD_ARGS=()
if [[ ${#METHODS[@]} -gt 0 ]]; then
    METHOD_ARGS=(--methods "${METHODS[@]}")
fi

echo "=== MolSimStack Run ==="
echo "Domain  : $DOMAIN"
echo "Input   : $INPUT"
echo "Output  : $OUTDIR"
echo "Threads : $THREADS"
echo "------------------------"

python "$BASE_DIR/core/dispatcher.py" \
    --input "$INPUT" \
    --outdir "$OUTDIR" \
    --domain "$DOMAIN" \
    --threads "$THREADS" \
    "${METHOD_ARGS[@]}"

echo "=== Completed ==="
