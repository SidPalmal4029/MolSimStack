import os
import numpy as np
from rdkit.Chem import Descriptors, rdMolDescriptors

from core import clustering, visualization
from core.utils import save_matrices

def compute_descriptor_vector(m):
    return np.array([
        Descriptors.MolWt(m),
        Descriptors.MolLogP(m),
        Descriptors.TPSA(m),
        Descriptors.NumHDonors(m),
        Descriptors.NumHAcceptors(m),
        Descriptors.NumRotatableBonds(m),
        Descriptors.HeavyAtomCount(m),
        rdMolDescriptors.CalcNumRings(m),
        rdMolDescriptors.CalcNumAromaticRings(m),
        rdMolDescriptors.CalcFractionCSP3(m)
    ])

#Standardization
def standardize(X):
    mean = np.mean(X, axis=0)
    std = np.std(X, axis=0)

    std[std == 0] = 1.0  # avoid division by zero

    return (X - mean) / std


def cosine_similarity_matrix(X):
    norm = np.linalg.norm(X, axis=1, keepdims=True)
    norm[norm == 0] = 1.0  # avoid division by zero
    Xn = X / norm
    return Xn @ Xn.T


def run(mols, names, outdir, threads=1, cluster_method="average"):
    os.makedirs(outdir, exist_ok=True)

    # Compute descriptor matrix
    X = np.array([compute_descriptor_vector(m) for m in mols])

    if len(mols) < 2:
        return

    # standardize features
    X_scaled = standardize(X)

    # Similarity
    sim = cosine_similarity_matrix(X_scaled)

    # Save + clustering
    dist = save_matrices(sim, names, outdir)
    np.fill_diagonal(dist, 0.0)

    link = clustering.hierarchical(dist)

    # Visualization
    visualization.plot_heatmap(
        sim, link, names,
        os.path.join(outdir, "heatmap.png"),
        method_name=f"Descriptors (Cosine, {cluster_method})"
    )
