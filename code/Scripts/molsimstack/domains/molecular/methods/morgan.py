import os
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

from ..similarity import tanimoto_matrix
from core import clustering, visualization
from core.utils import save_matrices

def run(mols, names, outdir, threads=1):
    os.makedirs(outdir, exist_ok=True)

    print("[Morgan] Generating fingerprints")

    gen = GetMorganGenerator(radius=2, fpSize=2048)
    fps = [gen.GetFingerprint(m) for m in mols]

    print("[Morgan] Computing similarity matrix")
    sim = tanimoto_matrix(fps)

    print("[Morgan] Clustering")
    dist = save_matrices(sim, names, outdir)
    link = clustering.hierarchical(dist)

    print("[Morgan] Plotting heatmap")
    visualization.plot_heatmap(
        sim,
        link,
        names,
        os.path.join(outdir, "heatmap.png"),
        method_name="Morgan (Tanimoto)"
    )

    print("[Morgan] Done")
