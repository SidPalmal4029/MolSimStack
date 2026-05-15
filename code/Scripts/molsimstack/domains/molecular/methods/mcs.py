import os
import numpy as np
import pandas as pd
from rdkit.Chem import rdFMCS

from core import clustering, visualization


def compute_pairwise_mcs_similarity(mols):
    n = len(mols)
    sim = np.zeros((n, n))

    for i in range(n):
        sim[i, i] = 1.0
        for j in range(i + 1, n):
            res = rdFMCS.FindMCS([mols[i], mols[j]], timeout=10)
            size = res.numAtoms

            size_i = mols[i].GetNumAtoms()
            size_j = mols[j].GetNumAtoms()

            score = size / min(size_i, size_j)

            sim[i, j] = score
            sim[j, i] = score

    return sim


def run(mols, names, outdir, threads=1):
    os.makedirs(outdir, exist_ok=True)

    print("[MCS] Computing global MCS")

    res = rdFMCS.FindMCS(
        mols,
        timeout=30,
        completeRingsOnly=True
    )

    smarts = res.smartsString

    with open(os.path.join(outdir, "mcs_global.txt"), "w") as f:
        f.write(smarts)

    print("[MCS] Global MCS saved")

    # Pairwise similarity
    print("[MCS] Computing pairwise similarity (slow)")

    sim = compute_pairwise_mcs_similarity(mols)

    df = pd.DataFrame(sim, index=names, columns=names)
    df.to_csv(os.path.join(outdir, "mcs_similarity.csv"))

    # Heatmap
    if len(mols) > 1:
        dist = 1 - sim
        np.fill_diagonal(dist, 0.0)

        link = clustering.hierarchical(dist)

        visualization.plot_heatmap(
            sim,
            link,
            names,
            outdir,
            method_key="mcs",
            method_name="MCS (Pairwise)"
        )

    print("[MCS] Done")
