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

    # -----------------------------
    # Find common molecules
    # -----------------------------
    common_mols = set(matrices[method_names[0]].index)
    for m in method_names[1:]:
        common_mols = common_mols.intersection(set(matrices[m].index))

    common_mols = sorted(common_mols)

    print(f"[PerMol] Processing {len(common_mols)} molecules across {len(method_names)} methods")

    results = []

    # -----------------------------
    # Iterate molecules
    # -----------------------------
    for mol in common_mols:
        for i in range(len(method_names)):
            for j in range(i + 1, len(method_names)):

                m1 = method_names[i]
                m2 = method_names[j]

                # safety check
                if mol not in matrices[m1].index or mol not in matrices[m2].index:
                    continue

                row1 = matrices[m1].loc[mol]
                row2 = matrices[m2].loc[mol]

                # -----------------------------
                # Align columns
                # -----------------------------
                common_cols = row1.index.intersection(row2.index)

                if len(common_cols) == 0:
                    continue

                # -----------------------------
                # Remove self safely
                # -----------------------------
                if mol in common_cols:
                    keep_cols = [c for c in common_cols if c != mol]
                else:
                    keep_cols = list(common_cols)

                if len(keep_cols) < 2:
                    continue

                v1 = row1[keep_cols].values
                v2 = row2[keep_cols].values

                # -----------------------------
                # Final shape safety
                # -----------------------------
                if len(v1) != len(v2):
                    continue

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

def summarize_per_molecule(df, outdir):
    comp_dir = os.path.join(outdir, "comparison")

    summary = df.groupby("Molecule")["Correlation"].mean().reset_index()
    summary.columns = ["Molecule", "MeanCorrelation"]

    out_file = os.path.join(comp_dir, "per_molecule_summary.csv")
    summary.to_csv(out_file, index=False)

    print(f"[PerMol] Summary saved to {out_file}")
