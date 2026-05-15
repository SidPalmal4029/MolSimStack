import os
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate

from ..similarity import tanimoto_matrix
from core import clustering, visualization
from core.utils import save_matrices

def run(mols, names, outdir, threads=1):
    os.makedirs(outdir, exist_ok=True)

    factory = Gobbi_Pharm2D.factory
    fps = [Generate.Gen2DFingerprint(m, factory) for m in mols]

    sim = tanimoto_matrix(fps)

    if len(mols) < 2:
        return

    dist = save_matrices(sim, names, outdir)
    link = clustering.hierarchical(dist)

    visualization.plot_heatmap(
        sim, link, names,outdir,
        method_key="pharmacophore",
        method_name="Pharmacophore"
    )
