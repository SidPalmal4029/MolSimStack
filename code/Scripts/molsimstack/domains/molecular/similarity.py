import numpy as np
from rdkit.DataStructs import TanimotoSimilarity

def tanimoto_matrix(fps):
    n = len(fps)
    mat = np.zeros((n, n), dtype=float)

    for i in range(n):
        mat[i, i] = 1.0
        for j in range(i + 1, n):
            sim = TanimotoSimilarity(fps[i], fps[j])
            mat[i, j] = sim
            mat[j, i] = sim

    return mat
