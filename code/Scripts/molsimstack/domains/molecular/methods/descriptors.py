import os
import numpy as np
from rdkit.Chem import Descriptors

from core import clustering, visualization
from core.utils import save_matrices

def compute_descriptor_vector(m):
    return np.array([
        Descriptors.MolWt(m),
        Descriptors.MolLogP(m),
        Descriptors.TPSA(m),
        Descriptors.NumHDonors(m),
        Descriptors.NumHAcceptors(m)
    ])


def cosine_similarity_matrix(X):
    norm = np.linalg.norm(X, axis=1, keepdims=True)
    Xn = X / norm
    return Xn @ Xn.T


def run(mols, names, outdir, threads=1):
    os.makedirs(outdir, exist_ok=True)

    X = np.array([compute_descriptor_vector(m) for m in mols])
    sim = cosine_similarity_matrix(X)

    if len(mols) < 2:
        return

    dist = save_matrices(sim, names, outdir)
    np.fill_diagonal(dist, 0.0)
    link = clustering.hierarchical(dist)

    visualization.plot_heatmap(
        sim, link, names,
        os.path.join(outdir, "heatmap.png"),
        method_name="Descriptors (Cosine)"
    )
