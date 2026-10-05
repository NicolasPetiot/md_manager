import MDAnalysis as mda
import networkx as nx
import numpy as np
from MDAnalysis.analysis import distances

from .core import universe_to_top


def get_inter_atomic_contacts(ag: mda.AtomGroup, cutoff: float, box=None) -> nx.Graph:
    G = nx.Graph()
    is_in_contact = distances.contact_matrix(ag.atoms.positions, cutoff=cutoff, box=box)
    edges = zip(*np.where(is_in_contact))  # pyright: ignore[reportArgumentType]
    G.add_edges_from(edges)

    node_attr = universe_to_top(ag).to_dict("index")
    nx.set_node_attributes(G, values=node_attr)

    return G
