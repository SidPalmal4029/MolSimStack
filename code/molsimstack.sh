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
    echo "MolSimStack - Multi-Method Molecular Similarity Pipeline"
    echo ""
    echo "Usage:"
    echo "  molsimstack.sh -i <input> -o <outdir> -d <domain> [options]"
    echo ""
    echo "Required:"
    echo "  -i, --input       Input file or directory"
    echo "                    mol → .sdf file OR directory of .sdf files"
    echo "                    nuc → FASTA (.fa / .fasta)"
    echo ""
    echo "  -o, --outdir      Output directory"
    echo ""
    echo "  -d, --domain      Domain:"
    echo "                    mol  → small molecules / peptides"
    echo "                    nuc  → nucleic acids (BLAST-based) future suport"
    echo "                    prot → proteins (future support)"
    echo ""
    echo "Optional:"
    echo "  -m, --methods     Methods to run (space-separated)"
    echo "                    Default (mol): morgan topo maccs shape3d"
    echo "                    Use 'all' to run all methods"
    echo ""
    echo "                    Available methods (mol):"
    echo "                      morgan         → circular fingerprints"
    echo "                      topo           → topological fingerprints"
    echo "                      maccs          → MACCS keys"
    echo "                      shape3d        → 3D shape similarity"
    echo "                      atompair       → atom pair fingerprints"
    echo "                      torsion        → torsion fingerprints"
    echo "                      descriptors    → physicochemical properties"
    echo "                      pharmacophore  → feature-based similarity"
    echo "                      scaffold       → Murcko scaffold grouping"
    echo "                      mcs            → maximum common substructure"
    echo ""
    echo "  -t, --threads     Number of threads (default: 1)"
    echo ""
    echo "Behavior:"
    echo "  • Multiple methods are executed sequentially."
    echo "  • Similarity and distance matrices are saved per method."
    echo "  • Clustering and heatmaps are generated where applicable."
    echo ""
    echo "Post-analysis (auto-enabled when ≥2 methods):"
    echo "  • Method comparison (global agreement)"
    echo "  • Per-molecule comparison (local consistency)"
    echo "  • Per-molecule summary (stability score)"
    echo ""
    echo ""
    echo "Examples:"
    echo "  # Default methods"
    echo "  molsimstack.sh -i molecules.sdf -o results -d mol"
    echo ""
    echo "  # Run selected methods"
    echo "  molsimstack.sh -i molecules.sdf -o results -d mol -m morgan topo maccs"
    echo ""
    echo "  # Run all methods"
    echo "  molsimstack.sh -i molecules.sdf -o results -d mol -m all"
    echo ""
    echo "  # Directory input"
    echo "  molsimstack.sh -i sdf_folder/ -o results -d mol"
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

LOGFILE="$OUTDIR/run.log"

mkdir -p "$OUTDIR"

echo "=== MolSimStack Run ===" | tee "$LOGFILE"
echo "Domain  : $DOMAIN"     | tee -a "$LOGFILE"
echo "Input   : $INPUT"      | tee -a "$LOGFILE"
echo "Output  : $OUTDIR"     | tee -a "$LOGFILE"
echo "Threads : $THREADS"    | tee -a "$LOGFILE"
echo "------------------------" | tee -a "$LOGFILE"

python "$BASE_DIR/core/dispatcher.py" \
    --input "$INPUT" \
    --outdir "$OUTDIR" \
    --domain "$DOMAIN" \
    --threads "$THREADS" \
    "${METHOD_ARGS[@]}" \
    2>&1 | tee -a "$LOGFILE"

echo "=== Completed ===" | tee -a "$LOGFILE"
