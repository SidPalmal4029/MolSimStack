import os
from rdkit.Chem import RDKFingerprint

from ..similarity import tanimoto_matrix
from core import clustering, visualization
from core.utils import save_matrices

def run(mols, names, outdir, threads=1):
    os.makedirs(outdir, exist_ok=True)

    print("[Topo] Generating fingerprints")

    fps = [RDKFingerprint(m) for m in mols]

    print("[Topo] Computing similarity matrix")
    sim = tanimoto_matrix(fps)

    if len(mols) < 2:
        print("[Topo] Only one molecule → skipping clustering")
        return

    print("[Topo] Clustering")
    dist = save_matrices(sim, names, outdir)
    link = clustering.hierarchical(dist)

    print("[Topo] Plotting heatmap")
    visualization.plot_heatmap(
        sim,
        link,
        names,
        outdir,
        method_key="topo",
        method_name="Topological (Tanimoto)"
    )

    print("[Topo] Done")
