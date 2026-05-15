import os
from rdkit.Chem.rdFingerprintGenerator import GetAtomPairGenerator

from ..similarity import tanimoto_matrix
from core import clustering, visualization
from core.utils import save_matrices

def run(mols, names, outdir, threads=1):
    os.makedirs(outdir, exist_ok=True)

    gen = GetAtomPairGenerator(fpSize=2048)
    fps = [gen.GetFingerprint(m) for m in mols]

    sim = tanimoto_matrix(fps)

    if len(mols) < 2:
        return

    dist = save_matrices(sim, names, outdir)
    link = clustering.hierarchical(dist)

    visualization.plot_heatmap(
        sim,
        link,
        names,
        outdir,
        method_key="atompair",
        method_name="AtomPair"
    )
