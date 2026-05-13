import os
import pandas as pd
import numpy as np
from rdkit.Chem.Scaffolds import MurckoScaffold

from core import clustering, visualization


def run(mols, names, outdir, threads=1):
    os.makedirs(outdir, exist_ok=True)

    # Generate scaffolds
    scaffolds = [MurckoScaffold.MurckoScaffoldSmiles(mol=m) for m in mols]

    # Save mapping
    df = pd.DataFrame({"Molecule": names, "Scaffold": scaffolds})
    df.to_csv(os.path.join(outdir, "scaffold_mapping.csv"), index=False)

    print("[Scaffold] Mapping saved")

    # Group by scaffold
    groups = df.groupby("Scaffold")["Molecule"].apply(list)

    with open(os.path.join(outdir, "scaffold_groups.txt"), "w") as f:
        for scaf, mol_list in groups.items():
            f.write(f"{scaf}\n")
            for m in mol_list:
                f.write(f"  - {m}\n")
            f.write("\n")

    print("[Scaffold] Groups saved")

    # Optional similarity matrix
    n = len(scaffolds)
    sim = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            sim[i, j] = 1.0 if scaffolds[i] == scaffolds[j] else 0.0

    df_sim = pd.DataFrame(sim, index=names, columns=names)
    df_sim.to_csv(os.path.join(outdir, "scaffold_similarity.csv"))

    print("[Scaffold] Similarity matrix saved")

    # Heatmap (optional)
    if n > 1:
        dist = 1 - sim
        link = clustering.hierarchical(dist)

        visualization.plot_heatmap(
            sim,
            link,
            names,
            os.path.join(outdir, "heatmap.png"),
            method_name="Scaffold (Binary)"
        )

    print("[Scaffold] Done")

