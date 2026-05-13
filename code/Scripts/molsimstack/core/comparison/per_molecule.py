import os
import pandas as pd
import numpy as np


def load_matrix(path):
    return pd.read_csv(path, index_col=0)


def per_molecule_comparison(outdir, methods):
    matrices = {}

    # -----------------------------
    # Load matrices
    # -----------------------------
    for m in methods:
        path = os.path.join(outdir, m, "similarity_matrix.csv")
        if os.path.exists(path):
            matrices[m] = load_matrix(path)
        else:
            print(f"[PerMol] Skipping {m} (no matrix)")

    if len(matrices) < 2:
        print("[PerMol] Not enough valid methods")
        return None

    method_names = list(matrices.keys())
    mols = matrices[method_names[0]].index

    print(f"[PerMol] Processing {len(mols)} molecules across {len(method_names)} methods")

    results = []

    # -----------------------------
    # Iterate molecules
    # -----------------------------
    for mol in mols:
        for i in range(len(method_names)):
            for j in range(i + 1, len(method_names)):

                m1 = method_names[i]
                m2 = method_names[j]

                v1 = matrices[m1].loc[mol].values
                v2 = matrices[m2].loc[mol].values

                # -----------------------------
                # Remove self index safely
                # -----------------------------
                idx = matrices[m1].index.get_loc(mol)

                v1 = np.delete(v1, idx)
                v2 = np.delete(v2, idx)

                # -----------------------------
                # Handle constant vectors
                # -----------------------------
                if np.std(v1) == 0 or np.std(v2) == 0:
                    corr = 0.0
                else:
                    corr = np.corrcoef(v1, v2)[0, 1]

                    if np.isnan(corr):
                        corr = 0.0

                results.append([mol, m1, m2, corr])

    df = pd.DataFrame(
        results,
        columns=["Molecule", "Method1", "Method2", "Correlation"]
    )

    return df


def save_per_molecule(df, outdir):
    comp_dir = os.path.join(outdir, "comparison")
    os.makedirs(comp_dir, exist_ok=True)

    out_file = os.path.join(comp_dir, "per_molecule_comparison.csv")
    df.to_csv(out_file, index=False)

    print(f"[PerMol] Saved to {out_file}")


# -----------------------------
# Optional: summary (very useful)
# -----------------------------
def summarize_per_molecule(df, outdir):
    comp_dir = os.path.join(outdir, "comparison")

    summary = df.groupby("Molecule")["Correlation"].mean().reset_index()
    summary.columns = ["Molecule", "MeanCorrelation"]

    out_file = os.path.join(comp_dir, "per_molecule_summary.csv")
    summary.to_csv(out_file, index=False)

    print(f"[PerMol] Summary saved to {out_file}")
