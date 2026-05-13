import os
from rdkit.Chem import MACCSkeys

from ..similarity import tanimoto_matrix
from core import clustering, visualization
from core.utils import save_matrices



def run(mols, names, outdir, threads=1):
    os.makedirs(outdir, exist_ok=True)

    print("[MACCS] Generating fingerprints")

    fps = [MACCSkeys.GenMACCSKeys(m) for m in mols]

    print("[MACCS] Computing similarity matrix")
    sim = tanimoto_matrix(fps)

    if len(mols) < 2:
        print("[MACCS] Only one molecule → skipping clustering")
        return

    print("[MACCS] Clustering")
    dist = save_matrices(sim, names, outdir)
    link = clustering.hierarchical(dist)

    print("[MACCS] Plotting heatmap")
    visualization.plot_heatmap(
        sim,
        link,
        names,
        os.path.join(outdir, "heatmap.png"),
        method_name="MACCS Keys (Tanimoto)"
    )

    print("[MACCS] Done")
