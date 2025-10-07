import networkx as nx
import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix

__all__ = [
    "frame_contact_graph",
    "frame_coarse_grained_contact_graph",
    "annotate_graph",
    "annotate_coarse_grained_graph",
    "get_contact_graph",
    "get_coarse_grained_contact_graph",
]

def frame_contact_graph(df:pd.DataFrame, threshold:float) -> nx.Graph:
    """
    
    """
    # Make sure indexes are range(Natm)
    Natm = len(df)
    idx = pd.Index(range(Natm), name="atom_id")
    df.index = idx

    xyz = ["x", "y", "z"]
    G = get_contact_graph(df[xyz], threshold=threshold)
    G = annotate_graph(G, top=df)
    return G

def frame_coarse_grained_contact_graph(df:pd.DataFrame, groups:list[str]|str, threshold:float) -> nx.Graph:
    """
    
    """

    # Make sure indexes are range(Natm)
    Natm = len(df)
    idx = pd.Index(range(Natm), name="atom_id")
    df.index = idx

    if isinstance(groups, str):
        groups = [groups]
    grouped_atoms = df.groupby(groups)

    # Build id2node & node2idx dictionaries:
    id2node  = {int(id): node for node, (_, df) in enumerate(grouped_atoms) for id in df.index}

    xyz = ["x", "y", "z"]
    G = get_coarse_grained_contact_graph(df[xyz], id2node=id2node, threshold=threshold)
    G = annotate_coarse_grained_graph(G, grouped_atoms, groups, top=df)
    return G


def get_contact_graph(pos1:np.ndarray|pd.DataFrame, pos2:np.ndarray|pd.DataFrame = None, threshold:float = 4.0) -> nx.Graph:
    """
    
    """
    is_inter = pos2 is not None # inter-structure contact graph (i.e. two distincts structures)

    # Check pos dimention
    _, Ndim = pos1.shape
    if Ndim != 3:
        raise ValueError("pos1 should be three dimentional")
    
    if is_inter:
        _, Ndim = pos2.shape
        if Ndim != 3:
            raise ValueError("pos2 should be three dimentional")
        
    else:
        # pos2 is None -> should be same as pos1
        pos2 = pos1.copy()        

    # pos1 and pos2 are three dimentional -> good for inter-node distance computations:
    dists = distance_matrix(pos1, pos2)

    # Build graph edges:
    I, J = np.where(dists < threshold)
    edges = list(zip(I, J))
    if not is_inter:
        # Should remove edges (i, i)
        edges = [(i, j) for i, j in edges if i != j]

    return nx.Graph(edges)

def annotate_graph(G:nx.Graph, top:pd.DataFrame) -> nx.Graph:
    """
    
    """
    Natm = len(top)
    if len(G) != Natm:
        raise ValueError(f"Graph contains {len(G)} nodes and Topology contains {Natm} atoms...")
    
    cols = top.columns

    # Atom Name:
    col = "name"
    if col in cols:
        d = {node: attr for node, attr in enumerate(top[col].values)}
        nx.set_node_attributes(G, d, col)

    # Residue index
    col = "resi"
    if col in cols:
        d = {node: attr for node, attr in enumerate(top[col].values)}
        nx.set_node_attributes(G, d, col)

    # Residue name
    col = "resn"
    if col in cols:
        d = {node: attr for node, attr in enumerate(top[col].values)}
        nx.set_node_attributes(G, d, col)

    # Chain labels:
    col = "chain"
    if col in cols:
        d = {node: attr for node, attr in enumerate(top[col].values)}
        nx.set_node_attributes(G, d, col)

        # Inter-chain contact: edge between two different chains
        d = {(i, j): d[i] != d[j] for i, j in G.edges}
        nx.set_edge_attributes(G, d, "is_inter_chain")

    return G

def annotate_coarse_grained_graph(G:nx.Graph, grouped_atoms, groups:list[str], top:pd.DataFrame):
    """
    
    """
    Ngrps = len(grouped_atoms)
    if len(G) != Ngrps:
        raise ValueError(f"Graph contains {len(G)} nodes and Topology contains {Ngrps} groups of atom...")
    
    group_labels = [cols for cols, _ in grouped_atoms]
    for i, attr in enumerate(groups):
        d = {node: label[i] for node, label in enumerate(group_labels)}
        nx.set_node_attributes(G, d, attr)

    return G


def get_coarse_grained_contact_graph(pos:np.ndarray|pd.DataFrame, id2node:dict[int, int], threshold:float = 4.0) -> nx.Graph:
    """
    
    """
    # Check pos dimention
    _, Ndim = pos.shape
    if Ndim != 3:
        raise ValueError("pos should be three dimentional")     

    dists = distance_matrix(pos, pos)

    # Build graph edges:
    I, J = np.where(dists < threshold)
    edges = list(zip(I, J))
    
    # Switch to coarse-grained representation:
    edges = [(id2node[i], id2node[j]) for i, j in edges]
    edges = [(i, j) for i, j in edges if i != j] # Should remove edges (i, i)
    return nx.Graph(edges)
