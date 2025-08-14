from .universe2df_relation import universe_to_top, top_to_universe

import MDAnalysis as mda
import pandas as pd

from typing import Self
from dataclasses import dataclass, asdict

__all__ = [
    "TrajParams",
    "Traj",
]

@dataclass
class TrajParams:
    # Topology Records:
    return_record_name:bool = True
    return_name:bool = True
    return_alt:bool = True
    return_resn:bool = True
    return_chain:bool = True
    return_resi:bool = True
    return_icode:bool = True
    return_occupancy:bool = True
    return_b:bool = True
    return_segi:bool = True
    return_e:bool = True
    return_q:bool = True
    return_m:bool = True
    return_type:bool = True
    return_atom_id:bool = True

    # Data Records:
    return_v:bool = False

class Traj:
    default_params = asdict(TrajParams())

    #################################################################################################################################################
    # Initialization code:
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
        self.universe = mda.Universe(*args, **kwargs)
        
        # Read Universe -> Top
        top = universe_to_top(self.universe)
        self.set_top(top)

    @classmethod
    def from_df(cls, df:pd.DataFrame, Nframe = 1, **kwargs) -> Self:
        """
        Initialization from DataFrame:

        Reads the topology and the coordinates stored in `df` and build the associated Traj object
        """
        # Params:
        traj = cls.__new__(cls)
        _ = traj.__init_params(**kwargs)

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

    #################################################################################################################################################
    # Allow slices:
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
        
    #################################################################################################################################################
    # Utils:
    
    def __len__(self):
        return len(self.universe.trajectory)

    def __repr__(self):
        return f"{self.__class__.__module__}.{self.__class__.__name__}: Nframe = {len(self)}; Natom = {len(self.top)}"
    
    def save(self, *topology, frames="all", **kwargs):
        self.universe.atoms.write(*topology, frames=frames, **kwargs)