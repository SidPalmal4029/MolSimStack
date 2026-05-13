import os
import pandas as pd
import numpy as np

def save_matrices(sim, names, outdir):
    os.makedirs(outdir, exist_ok=True)

    # -----------------------------
    # Similarity matrix
    # -----------------------------
    df_sim = pd.DataFrame(sim, index=names, columns=names)
    df_sim.to_csv(os.path.join(outdir, "similarity_matrix.csv"))

    # -----------------------------
    # Distance matrix
    # -----------------------------
    dist = 1 - sim
    np.fill_diagonal(dist, 0.0)

    df_dist = pd.DataFrame(dist, index=names, columns=names)
    df_dist.to_csv(os.path.join(outdir, "distance_matrix.csv"))

    return dist
