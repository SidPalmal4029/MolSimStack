from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

def hierarchical(distance_matrix):
    condensed = squareform(distance_matrix)
    return linkage(condensed, method="average")
