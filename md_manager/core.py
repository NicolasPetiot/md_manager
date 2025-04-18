from .parameters import TrajParams, ATTRIBUTE_RECORD_EQUIVALENCE

from MDAnalysis import Universe, NoDataError

from pathlib import Path
import pandas as pd

from dataclasses import asdict
from warnings import warn

import sys
from urllib.request import urlopen
from os.path import splitext

__all__ = ["Traj", "universe_to_top", "top_to_universe", "PDB", "GRO", "XYZ", "XTC", "load", "save", "fetch_PDB"]

def load(filename:str, *args, **kwargs) -> pd.DataFrame:    
    return Traj(filename, *args, **kwargs).load()

def save(df:pd.DataFrame, *filenames, frames = "all", **kwargs):
    Traj.from_df(df).save(*filenames, frames=frames, **kwargs)

class Traj:
    default_params = asdict(TrajParams())

    def __init__(self, *args, **kwargs):
        """
        Initialization code from file(s) .

        The system always requires a *topology file* --- in the simplest case just
        a list of atoms. This can be a CHARMM/NAMD PSF file or a simple coordinate
        file with atom informations such as XYZ, PDB, GROMACS GRO or TPR, or CHARMM
        CRD. See :ref:`Supported topology formats` for what kind of topologies can
        be read.

        A *trajectory file* provides coordinates; the coordinates have to be
        ordered in the same way as the list of atoms in the topology. A trajectory
        can be a single frame such as a PDB, CRD, or GRO file, or it can be a MD
        trajectory (in CHARMM/NAMD/LAMMPS DCD, GROMACS XTC/TRR, AMBER nc, generic
        XYZ format, ...).  See :ref:`Supported coordinate formats` for what can be
        read as a "trajectory".

        As a special case, when the topology is a file that contains atom
        information *and* coordinates (such as XYZ, PDB, GRO or CRD, see
        :ref:`Supported coordinate formats`) then the coordinates are immediately
        loaded from the "topology" file unless a trajectory is supplied.

        Parameters
        ----------
        """
        # Params:
        kwargs = self.__init_params(**kwargs)

        # Read File(s) -> Universe
        self.universe = Universe(*args, **kwargs)
        
        # Read Universe -> Top
        top = universe_to_top(self.universe)
        self.set_top(top)

    @classmethod
    def from_df(cls, df:pd.DataFrame, Nframe = 1, **kwargs):
        # Params:
        traj = cls.__new__(cls)
        __all__ = traj.__init_params(**kwargs)

        # Read df -> Top:
        traj.set_top(df)

        # Read Top -> Universe
        traj.universe = top_to_universe(traj.top, Nframe=Nframe)
        try:
            traj[0] = df

        except KeyError:
            pass

        return traj
    
    def __init_params(self, **kwargs) -> dict:
        """
        Use kwargs to update default params and set the parameters as attribute
        """
        params = self.default_params
        
        for key, val in kwargs.items():
            if key in params:
                params[key] = val

        # Set all attributes:
        for attr, val in params.items():
            setattr(self, attr, val)

        # Remove arguments that are already used:
        kwargs = {key: arg for key, arg in kwargs.items() if not key in params}

        return kwargs
    
    def set_top(self, top:pd.DataFrame) :
        """
        Use input DataFrame to create Traj Topology.
        """
        top_columns = [col.removeprefix("return_") for col, to_return in self.__dict__.items() if col.startswith("return_") and to_return]
        top_columns = [col for col in top_columns if col in top]
        self.top = top[top_columns]

    def __getitem__(self, idx) -> pd.DataFrame:
        """
        Returns the frame(s) indexed with idx.
        """
        if isinstance(idx, int):
            return self.load(idx)

        if isinstance(idx, slice):
            start = idx.start if idx.start is not None else 0
            stop  = idx.stop  if idx.stop  is not None else len(self)
            step  = idx.step  if idx.step  is not None else 1
            idx = list(range(start, stop, step))

        # check if idx is iterable:
        if hasattr(idx, '__iter__') and not isinstance(idx, str):
            return iter([self.load(id) for id in idx])
        
        raise ValueError("idx must be an int, slice or iterator")
    
    def __setitem__(self, id:int, df:pd.DataFrame):
        try:
            id = int(id)
            self.universe.trajectory[id].positions = df[["x", "y", "z"]].values
            if "vx" in df and "vy" in df and "vz" in df:
                self.universe.trajectory[id].velocities = df[["vx", "vy", "vz"]].values

        except TypeError:
            raise IndexError("Only support interger-like indices")

    
    def __len__(self):
        return len(self.universe.trajectory)

    def __repr__(self):
        return f"{self.__class__.__module__}.{self.__class__.__name__}: Nframe = {len(self)}; Natom = {len(self.top)}"
    
    def load(self, id:int = 0) -> pd.DataFrame:
        """
        Returns the DataFrame associated to the frame of index id
        """
        timestep = self.universe.trajectory[id] # Raises IndexError if id > length ot traj
        df = self.top.copy()

        # Position
        df[["x", "y", "z"]] = timestep.positions
        
        # Velocity:
        if self.return_v:
            df[["vx", "vy", "vz"]] = timestep.velocities
        
        return df
    
    def save(self, *topology, frames="all", **kwargs):
        self.universe.atoms.write(*topology, frames=frames, **kwargs)

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



class PDB(Traj):
    default_params = asdict(TrajParams(return_q = False, return_type = False, return_segi = False, return_icode = False))

class GRO(Traj):
    default_params = asdict(TrajParams(
        return_record_name = False,
        return_alt = False,
        return_chain = False,
        return_icode = False,
        return_occupancy = False,
        return_b = False,
        return_segi = False,
        return_q = False,
        return_type = False
    ))


class XYZ(Traj):
        default_params = asdict(TrajParams(
        # Topology Records:
        return_record_name = False,
        return_alt = False,
        return_resn = False,
        return_chain = False,
        return_resi = False,
        return_icode = False,
        return_occupancy = False,
        return_b = False,
        return_segi = False, 
        return_e = False,
        return_q = False,
        return_m = False,
        return_type = False,

        # Data Record:
        return_v = False,
    ))
        
class XTC(Traj):
    default_params = asdict(TrajParams(
        # Topology Records:
        return_record_name = False,
        return_alt = False,
        return_chain = False,
        return_icode = False,
        return_occupancy = False,
        return_b = False,
        return_segi = False, 
        return_e = False,
        return_m = False,
        return_type = False,

        # Data Record:
        return_v = False
    ))

def lines2df(lines:list[str], atom_only = False) -> pd.DataFrame:
    def read_format(line:str) -> tuple[str, str, str, str, int, str, float, float, float, float, float, float, str, str, str]:
        """
        Parses a line from a PDB file to extract atomic data.

        Args:
            line (str): A single line from the PDB file, formatted according 
                to the PDB specification.

        Returns:
            tuple[str, str, str, str, int, str, float, float, float, float, 
            float, float, str, str, str]: A tuple containing atom information, 
            including:
                - Atom name
                - Alternate location indicator
                - Residue name
                - Chain identifier
                - Residue sequence number
                - Insertion code
                - X, Y, and Z coordinates
                - Occupancy and temperature factor
                - Segment identifier, element symbol, and charge
        """
        return (
            line[:6].strip(),
            #int(line[6:11].strip()),
            line[12:16].strip(),        # Atom name
            line[16:17].strip(),        # Alternate location indicator
            line[17:20].strip() ,       # Residue name
            line[21:22].strip(),        # Chain identifier
            int(line[22:26].strip()),   # Residue sequence number
            line[26:27].strip(),        # Insertion code
            float(line[30:38].strip()), # X coordinate
            float(line[38:46].strip()), # Y coordinate
            float(line[46:54].strip()), # Z coordinate
            float(line[54:60].strip()) if line[54:60].strip() != "" else 1.0, # Occupancy
            float(line[60:66].strip()) if line[60:66].strip() != "" else 0.0, # Temperature factor
            line[72:76].strip(),        # Segment identifier
            line[76:78].strip(),        # Element symbol
            line[78:80].strip()         # Charge on the atom
        )
    
    columns = ["record_name", "name", "alt", "resn", "chain", "resi", "insertion", "x", "y", "z", "occupancy", "b", "segi", "e", "q"]
    data = [read_format(line) for line in lines if line[:6].strip() == "ATOM" or not atom_only and line[:6].strip() == "HETATM"]
    df = pd.DataFrame(data)
    df.columns = columns
    df.index.name = "atom_id"
    df.index += 1
    return df


def fetch_PDB(pdb_code:str, atom_only = False) -> pd.DataFrame:
    """
    Fetches a PDB structure from the RCSB PDB website and returns it as a pandas DataFrame.

    This function downloads a PDB file from the RCSB website based on the provided PDB code, 
    parses the structure, and returns it as a DataFrame. The DataFrame includes atomic data 
    from the PDB file. Optionally, only ATOM records can be returned.

    Args:
        pdb_code (str): The PDB code of the structure to fetch. This is a four-character identifier
                         assigned to each PDB entry.
        atom_only (bool, optional): If `True`, the returned DataFrame will only include 'ATOM' records
                                    (ignoring any HETATM or other record types). Defaults to `False`.

    Returns:
        pd.DataFrame: A DataFrame containing the atomic information from the fetched PDB file. 
                      The columns correspond to the fields defined by the `PDB` class.

    Raises:
        HTTPError: If the PDB code is invalid or the RCSB PDB website is unreachable.
        ValueError: If the PDB file cannot be parsed properly.

    Example:
        >>> df = fetch_PDB("1A2B", atom_only=True)
        >>> print(df.head())
           record_name name alt resn chain  resi insertion      x      y      z occupancy    b   segi e  q
        1          ATOM  CA   ALA     A    23          A  1.234  2.345  3.456     1.00  20.5  A  C  0
        2          ATOM  C    ALA     A    23          A  2.234  3.345  4.456     1.00  18.5  A  C  0
    """
    url = f"https://files.rcsb.org/download/{pdb_code.lower()}.pdb"
    response = urlopen(url)

    txt = response.read()
    lines = (txt.decode("utf-8") if sys.version_info[0] >= 3 else txt.decode("ascii")).splitlines()

    df = lines2df(lines)

    if atom_only:
        df = df.query("record_name == 'ATOM'")

    return df