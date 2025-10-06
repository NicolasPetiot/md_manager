from .parameters import ATTRIBUTE_RECORD_EQUIVALENCE

from MDAnalysis import Universe, NoDataError

import pandas as pd

from warnings import warn

__all__ = ["universe_to_top", "top_to_universe"]

def universe_to_top(u:Universe) -> pd.DataFrame:
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

def top_to_universe(top:pd.DataFrame, Nframe = 1) -> Universe:
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


