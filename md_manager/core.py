from MDAnalysis import Universe, NoDataError
from MDAnalysis.coordinates.PDB import PDBWriter
from MDAnalysis.coordinates.GRO import GROWriter
from MDAnalysis.coordinates.XYZ import XYZWriter
from MDAnalysis.coordinates.XTC import XTCWriter

from pathlib import Path
import pandas as pd

from os.path import splitext

__all__ = ["Traj", "PDB", "GRO", "XYZ", "XTC", "load", "save"]

def load(filename:str, *args, **kwargs) -> pd.DataFrame:
    file_type = {
        ".pdb":PDB,
        ".gro":GRO,
        ".xyz":XYZ,
    }
    _, ext = splitext(filename)
    if ext in file_type:
        return file_type[ext](filename, *args, **kwargs).load()
    
    if len(args) > 0 and splitext(args[0][1]) == ".xtc":
        return XTC(filename, *args, **kwargs).load()
    
    return Traj(filename, *args, **kwargs).load()

def save(df:pd.DataFrame, *filenames, **kwargs):
    file_type = {
        ".pdb":PDB,
        ".gro":GRO,
        ".xyz":XYZ,
        ".xtc":XTC
    }

    ext = splitext(filenames[-1])[1]
    try :
        traj = file_type[ext](df, **kwargs)
        traj.save(*filenames)

    except KeyError:
        raise ValueError(f"Unknown format {ext}. Must be in {list(file_type)}")


class Traj:
    default_params = dict(
        # Topology Records:
        return_record_name = False,
        return_name = True,
        return_alt = False,
        return_resn = True,
        return_chain = False,
        return_resi = True,
        return_icode = False,
        return_occupancy = False,
        return_b = False,
        return_segi = False, 
        return_e = False,
        return_q = False,
        return_m = False,
        return_type = False,
        return_atom_id = True,

        # Data Record:
        return_v = False,

        writer = None,
    )
    __attribute_record_equivalence = [
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

    # Dunder methods:
    # __init__
    # __getitem__
    # __setitem__
    # __len__
    # __str__

    def __init__(self, topology, *args, **kwargs):
        self.col2attr = {col: attr for attr, col in self.__attribute_record_equivalence}
        self.attr2col = {attr: col for attr, col in self.__attribute_record_equivalence}

        # Use kwargs to initialize the parameters
        kwargs = self.__init_params(**kwargs)

        # Initialize Traj from topology:
        if isinstance(topology, pd.DataFrame):
            self.__init_from_df(topology, **kwargs)
            return 
        
        try:
            topology = Path(topology)
            if topology.is_file() and topology.exists():
                self.__init_from_path(topology, *args, **kwargs)
                return

            else:
                raise FileNotFoundError(f"{topology} does not exists or is not acessible...")

        except (ValueError, TypeError):
            raise TypeError(f"topology should PathLike or a DataFrame...")
        
    def __getitem__(self, idx) -> pd.DataFrame:
        if isinstance(idx, int):
            return self.load(idx)

        if isinstance(idx, slice):
            start = idx.start if idx.start is not None else 0
            stop  = idx.stop  if idx.stop  is not None else len(self)
            step  = idx.step  if idx.step  is not None else 1
            idx = list(range(start, stop, step))

        # check if idx is iterable:
        if hasattr(idx, '__iter__') and not isinstance(idx, str):
            df =  pd.concat([self.load(i) for i in idx], ignore_index=True)
            df.index = pd.MultiIndex.from_product([idx, self.top.index], names=["frame", "atom_id"])
            return df


        raise TypeError("Index must be an integer, slice, or iterable of integers.")
    
    def __setitem__(self, id, df):
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
    
    # Utils:
    # load
    # save
    # close

    def load(self, id:int = 0) -> pd.DataFrame:
        """
        Returns the dataframe associated to the frame of index id
        """
        timestep = self.universe.trajectory[id] # Raises IndexError if id > length ot traj
        df = self.top.copy()

        # Position
        df[["x", "y", "z"]] = timestep.positions
        
        # Velocity:
        if self.return_v:
            df[["vx", "vy", "vz"]] = timestep.velocities
        
        return df
    
    def save(self, *filename, **kwargs):
        if self.writer is None:
            raise ValueError("No specified writer...")
        
        with self.writer(*filename, **kwargs) as file:
            file.write(self.universe)

    def close(self):
        self.universe.trajectory.close()

    # Various initialization methods:
    # __init_params 
    # __init_from_path
    # __init_from_df
    # __set_top

    def __init_params(self, **kwargs) -> dict:
        """
        Use kwargs to update default params and set the parameters as attribute
        """
        params = self.default_params.copy()
        self.columns = []
        for key, val in kwargs.items():
            if key in params:
                params[key] = val

        # Set all attributes:
        for attr, val in params.items():
            setattr(self, attr, val)
        
        self.columns    = [col for _, col in self.__attribute_record_equivalence if getattr(self, "return_" + col)]
        self.attributes = [self.col2attr[col] for col in self.columns]

        # Remove arguments that are already used:
        kwargs = {key: arg for key, arg in kwargs.items() if not key in params}

        return kwargs
    
    def __init_from_path(self, topology:Path, *args, **kwargs) -> None:
        self.universe = Universe(topology, *args, **kwargs)

        self.top = {}
        for attr in self.attr2col:
            try:
                self.top[attr] = getattr(self.universe.atoms, attr)

            except NoDataError:
                pass
        self.top = pd.DataFrame(self.top)
        self.top.columns = [self.attr2col[attr] for attr in self.top.columns]
        self.top = self.top[[col for col in self.top.columns if getattr(self, f"return_{col}")]]

        if "atom_id" in self.top:
            self.top = self.top.set_index("atom_id")

    def __init_from_df(self, df:pd.DataFrame, **kwargs):
        self.__set_top(df, Nframe=kwargs["Nframe"] if "Nframe" in kwargs else 1)
        self.__setitem__(id=0, df=df)

    def __set_top(self, top:pd.DataFrame, Nframe = 1):
        top = top[[col for col in top if col in self.col2attr]]
        residue_groups = [col for col in ["record_name", "chain", "resi"] if col in top]

        Natm = len(top)
        Nres = len(top.groupby(residue_groups)) if len(residue_groups) > 0 else 1
        atom_residx = top.resi.values - 1 if "resi" in top else top.index.values

        if Nres > 1:
            self.universe = Universe.empty(n_atoms=Natm, n_residues=Nres, atom_resindex=atom_residx, n_frames=Nframe, trajectory=True)

        else:
            self.universe = Universe.empty(n_atoms=Natm, n_frames=Nframe, trajectory=True)

        # Topology Attributes:
        for col in top:
            attr = self.col2attr[col]
            if not attr in {"resnames", "resids", "icodes", "segindices"}:
                self.universe.add_TopologyAttr(attr, top[col].values)

            elif attr == "resnames":
                residue_groups = [col for col in ["record_name", "chain", "resi"] if col in top]
                resn = top.groupby(residue_groups).resn.apply(lambda s: s.unique()[0]).values
                self.universe.add_TopologyAttr("resnames", resn)

            elif attr == "resids":
                resi = top.groupby(residue_groups).resi.apply(lambda s: s.unique()[0]).values
                self.universe.add_TopologyAttr("resids", resi)
        self.top = top

class PDB(Traj):
    default_params = dict(
        # Topology Records:
        return_record_name = True,
        return_name = True,
        return_alt = False,
        return_resn = True,
        return_chain = True,
        return_resi = True,
        return_icode = False,
        return_occupancy = False,
        return_b = True,
        return_segi = True, 
        return_e = False,
        return_q = False,
        return_m = False,
        return_type = False,
        return_atom_id = True,

        # Data Record:
        return_v = False,

        writer = PDBWriter,
    )

class GRO(Traj):
    default_params = dict(
        # Topology Records:
        return_record_name = False,
        return_name = True,
        return_alt = False,
        return_resn = True,
        return_chain = False,
        return_resi = True,
        return_icode = False,
        return_occupancy = False,
        return_b = False,
        return_segi = False, 
        return_e = False,
        return_q = False,
        return_m = False,
        return_type = False,
        return_atom_id = True,

        # Data Record:
        return_v = False,

        writer = GROWriter,
    )

class XYZ(Traj):
        default_params = dict(
        # Topology Records:
        return_record_name = False,
        return_name = True,
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
        return_atom_id = True,

        # Data Record:
        return_v = False,

        writer = XYZWriter,
    )
        
class XTC(Traj):
    default_params = dict(
        # Topology Records:
        return_record_name = False,
        return_name = True,
        return_alt = False,
        return_resn = True,
        return_chain = False,
        return_resi = True,
        return_icode = False,
        return_occupancy = False,
        return_b = False,
        return_segi = False, 
        return_e = False,
        return_q = True,
        return_m = False,
        return_type = False,
        return_atom_id = True,

        # Data Record:
        return_v = False,

        writer = XTCWriter,
    )