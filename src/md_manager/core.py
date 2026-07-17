import MDAnalysis as mda
import pandas as pd

from MDAnalysis import NoDataError

__all__ = ["universe_to_top"]

ATTRIBUTE_RECORD_EQUIVALENCE = [
    ("record_types", "record_name"),
    ("names", "name"),
    ("altLocs", "alt"),
    ("resnames", "resn"),
    ("chainIDs", "chain"),
    ("resids", "resi"),
    ("icodes", "icode"),
    ("occupancies", "occupancy"),
    ("tempfactors", "b"),
    ("segindices", "segi"),
    ("elements", "e"),
    ("charges", "q"),
    ("masses", "m"),
    ("types", "type"),
    ("ids", "atom_id")
]

class Traj(mda.Universe):
    def __init__(self, topology=None, *coordinates, all_coordinates=False, format=None, topology_format=None, transformations=None, guess_bonds=False, vdwradii=None, fudge_factor=0.55, lower_bound=0.1, in_memory=False, context="default", to_guess=[], force_guess=[], in_memory_step=1, **kwargs):
        super().__init__(topology, *coordinates, all_coordinates=all_coordinates, format=format, topology_format=topology_format, transformations=transformations, guess_bonds=guess_bonds, vdwradii=vdwradii, fudge_factor=fudge_factor, lower_bound=lower_bound, in_memory=in_memory, context=context, to_guess=to_guess, force_guess=force_guess, in_memory_step=in_memory_step, **kwargs)
        self.top = universe_to_top(self.universe)

def universe_to_top(u:mda.Universe) -> pd.DataFrame:
    """
    Function used to read the topology attributes of an input Universe and create the associated DataFrame.

    For all attributes in the ATTRIBUTE_RECORD_EQUIVALENCE list, the function tries to get the data from the Universe object and convert it to a Pandas Series.
    """
    top = {}
    for attr, col in ATTRIBUTE_RECORD_EQUIVALENCE:
        try:
            top[col] = getattr(u.atoms, attr)

        except NoDataError:
            pass

    top = pd.DataFrame(top)
    if "atom_id" in top:
        top = top.set_index("atom_id")
    return top
