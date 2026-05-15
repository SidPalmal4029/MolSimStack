import os
import numpy as np
from rdkit.Chem import AllChem
from rdkit.Chem import rdShapeHelpers

from core import clustering, visualization
from core.utils import save_matrices

def generate_3d(mol):
    mol = AllChem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)
    return mol


def shape_similarity_matrix(mols):
    n = len(mols)
    mat = np.zeros((n, n))

    for i in range(n):
        mat[i, i] = 1.0
        for j in range(i + 1, n):
            sim = 1 - rdShapeHelpers.ShapeTanimotoDist(mols[i], mols[j])
            mat[i, j] = sim
            mat[j, i] = sim

    return mat


def run(mols, names, outdir, threads=1):
    os.makedirs(outdir, exist_ok=True)

    print("[Shape3D] Generating 3D conformers")

    mols_3d = [generate_3d(m) for m in mols]

    print("[Shape3D] Computing similarity matrix")
    sim = shape_similarity_matrix(mols_3d)

    if len(mols) < 2:
        print("[Shape3D] Only one molecule → skipping clustering")
        return

    print("[Shape3D] Clustering")
    dist = save_matrices(sim, names, outdir)
    link = clustering.hierarchical(dist)

    print("[Shape3D] Plotting heatmap")
    visualization.plot_heatmap(
        sim,
        link,
        names,
        outdir,
        method_key="shape3d",
        method_name="3D Shape (Tanimoto)"
    )

    print("[Shape3D] Done")
