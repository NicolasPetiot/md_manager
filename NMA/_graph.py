from .._parameters import *
import networkx as nx
import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix
from scipy.linalg import eigh

### Create Graph Representation ###
def df2graph(df:pd.DataFrame, contact_threshold = 4.0) -> nx.Graph:
    """
    Returns a networkx.Graph object that contains nodes associated to all the atoms in the input df and that are connected if the distances between the atoms is smaller than an input contact_threshold.
    """
    G = nx.Graph()
    xyz = df[["x", "y", "z"]].to_numpy(dtype=float)
    distance_inter_nodes = distance_matrix(xyz, xyz)
    I, J = np.where((distance_inter_nodes < contact_threshold) & (distance_inter_nodes > 0.0))
    pairs = [(i, j) for i, j in zip(I, J)]
    G.add_edges_from(pairs)
    return G

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