import seaborn as sns
import matplotlib.pyplot as plt

def plot_heatmap(sim_matrix, linkage_matrix, labels, outpath, method_name=""):
    g = sns.clustermap(
        sim_matrix,
        row_linkage=linkage_matrix,
        col_linkage=linkage_matrix,
        xticklabels=labels,
        yticklabels=labels,
        cmap="viridis",
        annot=True,              # ✅ show values
        fmt=".2f",               # ✅ 2 decimal places
        linewidths=0.5
    )

    # -----------------------------
    # Add title (method label)
    # -----------------------------
    title = f"Similarity Heatmap ({method_name})" if method_name else "Similarity Heatmap"
    g.fig.suptitle(title, y=1.05)

    plt.savefig(outpath, bbox_inches="tight")
    plt.close()
