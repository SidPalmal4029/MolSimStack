#!/usr/bin/env bash
set -e

ENV_NAME="molsimstack"
PYTHON_VERSION="3.12"

echo "=== MolSimStack Installer ==="
# 1. Check conda
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found. Install Miniconda/Anaconda first."
    exit 1
fi

# Initialize conda for non-interactive shells
source "$(conda info --base)/etc/profile.d/conda.sh"

# 2. Detect mamba (or install it)
SOLVER="conda"

if command -v mamba &> /dev/null; then
    SOLVER="mamba"
    echo "Using mamba for faster installs"
else
    echo "mamba not found, installing into base..."
    conda install -y -n base -c conda-forge mamba
    SOLVER="mamba"
fi

# 3. Create environment
echo "Creating environment: $ENV_NAME"

$SOLVER create -y -n $ENV_NAME \
    python=$PYTHON_VERSION \
    numpy scipy pandas matplotlib seaborn scikit-learn \
    -c conda-forge


# 4. Activate environment
conda activate $ENV_NAME

# 5. Install RDKit (latest)
echo "Installing latest RDKit..."

$SOLVER install -y rdkit -c conda-forge

# 6. Install optional Python tools
echo "Installing CLI and utility packages..."

pip install typer pydantic loguru joblib

# 7. Verify installation
echo "Verifying RDKit..."

python - <<EOF
from rdkit import Chem
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

mol = Chem.MolFromSmiles("CCO")
gen = GetMorganGenerator(radius=2, fpSize=2048)
fp = gen.GetFingerprint(mol)

print("RDKit OK | FP length:", len(fp))
EOF

echo "=== Installation Complete ==="
echo "Activate environment with:"
echo "conda activate $ENV_NAME"

# -----------------------------
# Install engine code
# -----------------------------
echo "[3/5] Installing core engine"

INSTALL_DIR="$CONDA_PREFIX/share/molsimstack"
mkdir -p "$INSTALL_DIR"

# Clean previous install
rm -rf "$INSTALL_DIR"/*

# Copy working code ONLY
cp -r code/Scripts/molsimstack/* "$INSTALL_DIR/"

# -----------------------------
# Install wrapper (CLI)
# -----------------------------
echo "[4/5] Installing CLI wrapper"

WRAPPER_SRC="code/molsimstack.sh"
BIN_PATH="$CONDA_PREFIX/bin/molsimstack"

# Check wrapper exists
if [[ ! -f "$WRAPPER_SRC" ]]; then
    echo "ERROR: CLI wrapper not found: $WRAPPER_SRC"
    echo "Please ensure molsimstack.sh exists in the repository."
    exit 1
fi

# Copy wrapper
cp "$WRAPPER_SRC" "$BIN_PATH"

# Make executable
chmod +x "$BIN_PATH"

# -----------------------------
# Done
# -----------------------------
echo "[5/5] Installation complete"

echo ""
echo "Activate:"
echo "conda activate $ENV_NAME"
echo ""
echo "Test:"
echo "molsimstack -h"
