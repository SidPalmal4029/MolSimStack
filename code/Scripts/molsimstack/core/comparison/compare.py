import os
import pandas as pd
import numpy as np


def load_matrix(path):
    return pd.read_csv(path, index_col=0).values


def compare_methods(outdir, methods):
    matrices = {}

    # -----------------------------
    # Load matrices
    # -----------------------------
    for m in methods:
        path = os.path.join(outdir, m, "similarity_matrix.csv")

        if os.path.exists(path):
            matrices[m] = load_matrix(path)
        else:
            print(f"[Comparison] Skipping {m} (no matrix)")

    if len(matrices) < 2:
        print("[Comparison] Less than 2 valid methods → skipping comparison")
        return None

    # -----------------------------
    # Compute differences
    # -----------------------------
    results = []

    method_list = list(matrices.keys())

    for i in range(len(method_list)):
        for j in range(i + 1, len(method_list)):
            m1 = method_list[i]
            m2 = method_list[j]

            diff = np.abs(matrices[m1] - matrices[m2])
            score = diff.mean()

            results.append((m1, m2, score))

    return results


def save_results(results, outdir):
    comp_dir = os.path.join(outdir, "comparison")
    os.makedirs(comp_dir, exist_ok=True)

    df = pd.DataFrame(results, columns=["Method1", "Method2", "MeanDifference"])

    out_file = os.path.join(comp_dir, "method_comparison.csv")
    df.to_csv(out_file, index=False)

    print(f"[Comparison] Results saved to {out_file}")
