import MDAnalysis as mda
import pandas as pd

from MDAnalysis import NoDataError
from warnings import warn

__all__ = ["Traj", "universe_to_top", "top_to_universe"]

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

    def __getitem__(self, i:int) -> pd.DataFrame:
        df = self.top.copy()
        self.trajectory[i]
        if hasattr(self.atoms, "positions"):
            df[["x", "y", "z"]] = getattr(self.atoms, "positions")

        if hasattr(self.atoms, "velocities"):
            df[["vx", "vy", "vz"]] = getattr(self.atoms, "velocities")

        if hasattr(self.atoms, "forces"):
            df[["fx", "fy", "fz"]] = getattr(self.atoms, "forces")

        return df

    def __len__(self) -> int:
        return len(self.trajectory)


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

def top_to_universe(top:pd.DataFrame, Nframe = 1) -> mda.Universe:
    """
    Function used to read the topology records of an input DataFrame and create the associated Universe.
    """
    Natm = len(top)

    # Residues:
    groups = ["record_name", "chain", "resi"]
    groups = [group for group in groups if group in top]
    if len(groups) == 0:
        Nres = 1
        residues = None
        atm_resindex = None

    else:
        residues = top.groupby(groups)
        Nres = len(residues)
        Natm_residues = residues.name.count().values
        atm_resindex = [i for i, n in enumerate(Natm_residues) for _ in range(n)]

    # Segments:
    groups = ["record_name", "chain"]
    groups = [group for group in groups if group in top]
    if len(groups) == 0:
        Nseg = 1
        segments = None
        res_segindex = None

    else:
        segments = top.groupby(groups)
        Nseg = len(segments)
        res_segindex=[i for i, (_, grp) in enumerate(segments) for _ in range(len(grp.resi.unique()))]

    u = Universe.empty(n_atoms=Natm, n_residues=Nres, n_segments=Nseg, n_frames=Nframe, atom_resindex=atm_resindex, residue_segindex=res_segindex, trajectory=True)

    for attr, col in ATTRIBUTE_RECORD_EQUIVALENCE:
        if col in top:
            if not attr in {"resnames", "resids", "icodes", "segindices"}:
                u.add_TopologyAttr(attr, top[col].values)

            elif attr == "resnames" and residues is not None:
                resn = residues.resn.apply(lambda s: s.unique()[0]).values
                u.add_TopologyAttr("resnames", resn)

            elif attr == "resids" and residues is not None:
                resi = residues.resi.apply(lambda s: s.unique()[0]).values
                u.add_TopologyAttr("resids", resi)

            else:
                warn(f"{col} records are not yet supported by `top_to_universe`")

    return u
