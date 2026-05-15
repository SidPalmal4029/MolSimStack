import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import to_tree


# -----------------------------
# Convert linkage → Newick
# -----------------------------
def linkage_to_newick(linkage, labels):
    tree = to_tree(linkage, rd=False)

    def build(node):
        if node.is_leaf():
            return labels[node.id]
        left = build(node.left)
        right = build(node.right)
        return f"({left},{right})"

    return build(tree) + ";"


# -----------------------------
# Main plotting function
# -----------------------------
def plot_heatmap(sim, link, names, outdir, method_key, method_name="", figsize=None):
    n = len(names)

    # -----------------------------
    # Dynamic figure scaling
    # -----------------------------
    if n <= 100:
        scale = 0.5
    else:
        scale = 0.35

    fig_w = min(max(12, n * scale), 50)
    fig_h = min(max(10, n * scale), 50)

    if figsize is not None:
        fig_w, fig_h = figsize

    # -----------------------------
    # Annotation + label rules
    # -----------------------------
    if n <= 100:
        annot = True
        fmt = ".2f"
        label_size = max(4, 12 - n // 10)
    else:
        annot = False
        fmt = ""
        label_size = max(2, 8 - n // 100)

    # -----------------------------
    # Plot heatmap
    # -----------------------------
    g = sns.clustermap(
        sim,
        row_linkage=link,
        col_linkage=link,
        cmap="viridis",
        figsize=(fig_w, fig_h),
        xticklabels=names,
        yticklabels=names,
        dendrogram_ratio=(0.12, 0.12),
        cbar_kws={"label": "Similarity"},
        annot=annot,
        fmt=fmt
    )

    # -----------------------------
    # Label formatting
    # -----------------------------
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=label_size)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=label_size)

    # -----------------------------
    # Title
    # -----------------------------
    plt.title(method_name, fontsize=12)

    # -----------------------------
    # Output paths
    # -----------------------------
    heatmap_path = os.path.join(outdir, f"{method_key}_heatmap.png")
    nwk_path = os.path.join(outdir, f"{method_key}_dendogram.nwk")

    # -----------------------------
    # Save heatmap
    # -----------------------------
    plt.savefig(heatmap_path, dpi=300)
    plt.close()

    # -----------------------------
    # Save dendrogram (Newick)
    # -----------------------------
    newick_str = linkage_to_newick(link, names)
    with open(nwk_path, "w") as f:
        f.write(newick_str)
