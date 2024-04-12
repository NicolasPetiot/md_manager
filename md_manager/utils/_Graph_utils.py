from ..utils._params import *
import networkx as nx
import numpy as np
from scipy.linalg import eigh

__all__ = ["graph_laplacian", "graph_collective_modes"]

### Graph Laplacian ###
def graph_laplacian(G:nx.Graph) -> np.ndarray:
    """
    Returns a np.ndarray that corresponds to the inputed graph's laplacian.
    """
    Nnode = len(G)
    L     = np.zeros((Nnode, Nnode))
    for i, j in G.edges :
        Lij = -1
        L[i, j] = Lij
        L[j, i] = Lij
        L[i, i] -= Lij
        L[j, j] -= Lij

    return L

def graph_collective_modes(laplacian):
    """
    Returns a tuple containing the eigenvalues and eigenfrequencies of the inputed laplacian. Note that the first null eigenmodes have been removed.
    """
    eigenvals, eigenvecs = eigh(laplacian)
    eigenvals = eigenvals[1:]
    eigenvecs = eigenvecs[:, 1:]
    
    return eigenvals, eigenvecs