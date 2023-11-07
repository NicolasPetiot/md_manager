from ._parameters import *
from ._AtomSelection import AtomSelection
import networkx as nx
import numpy as np
from scipy.spatial import distance_matrix
from scipy.linalg import eigh

### Create Graph Representation ###
def sele2graph(sele:AtomSelection, contact_threshold = 4.0):
    """
    
    """
    nodes_position = sele.get_xyz()
    distance_inter_nodes = distance_matrix(nodes_position, nodes_position)
    Nnode = len(sele)

    G = nx.Graph()
    for i in range(Nnode):
        G.add_node(i)
        for j in range(i):
            if i != j and distance_inter_nodes[i, j] <= contact_threshold:
                G.add_edge(i, j)
    return G

### Graph Laplacian ###
def graph_laplacian(G:nx.Graph):
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