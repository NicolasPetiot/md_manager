from ..parameters import ONE_LETTER_CODE

from MDAnalysis import Universe
import pandas as pd

class Traj:
    topology_records = ["record_name", "name", "resn", "resi"]
    coord_records    = [ "x",  "y",  "z"]
    velocity_records = ["vx", "vy", "vz"]
    #force_records    = ["fx", "fy", "fz"]

    default_params = dict(
        # Read/Write mode:
        mode = "r",

        # Topology records:
        return_topology = True,     
        return_alt = False,
        return_chain = False,
        return_occupancy = False,
        return_b = False,
        return_segi = False,
        return_e = False,
        return_q = False,
        return_m = False,

        # Data records
        return_velocity = False,
        #return_forces = False # TODO

        # Secondary structures:
        return_alpha = False,
        return_beta = False,
        ss_method = "CUTABI"
    )

    def __init__(self, *args, **kwargs):
        # Set all default parameters with modifications according to kwargs
        self.reset(**kwargs)
        kwargs = {key: param for key, param in kwargs.items() if not key in self.default_params}

        # Initialize Universe Object: 
        if self.mode == "r":
            self.universe = Universe(*args, **kwargs)

        else:
            df = kwargs["df"] if "df" in kwargs else args[-1]
            if not isinstance(df, pd.DataFrame):
                raise ValueError("In write mode, please specify a pandas DataFrame...")            
            self.universe = self.__universe_from_df__(df)

        
        self.atoms = self.universe.atoms
        if self.return_topology:
            self.top = self.get_topology_df()
        

    def __universe_from_df__(self, df:pd.DataFrame) -> Universe:
        """
        
        """
        ca = df.query("name == 'CA'")
        Natm = len(df)
        Nres = len(ca)

        # Fill empy Universe with topology infos
        u = Universe.empty(Natm, Nres, atom_resindex=df.resi.values, trajectory=True)
        u.add_TopologyAttr("names", df.name.values)

        if "resi" in self.topology_records:
            u.add_TopologyAttr("resids", ca.resi.values)

        if "resn" in self.topology_records:
            u.add_TopologyAttr("resnames", ca.resn.values)

        if self.return_alt:
            u.add_TopologyAttr("altLocs", df["alt"].values)

        if self.return_chain:
            u.add_TopologyAttr("chainIDs", df["chain"].values)

        if self.return_occupancy:
            u.add_TopologyAttr("occupancies", df["occupancy"].values)

        if self.return_b:
            u.add_TopologyAttr("tempFactors", df["b"].values)

        if self.return_segi:
            u.add_TopologyAttr("segids", df["segi"].unique())

        if self.return_e:
            u.add_TopologyAttr("elements", df["e"].values)

        if self.return_q:
            u.add_TopologyAttr("charges", df["q"].values)

        if self.return_m:
            u.add_TopologyAttr("masses", df["m"].values)

        u.atoms.positions = df[self.coord_records].values

        return u

    def __len__(self):
        return len(self.universe.trajectory)        

    def __iter__(self):
        self.frame = iter(self.universe.trajectory)
        return self

    def __next__(self):
        u = next(self.frame)
        if self.return_topology:
            df = self.top.copy()
        else:
            df = pd.DataFrame(index=pd.Index(range(1, len(u)+1), name="atom_id"))
        

        df[self.coord_records] = u.positions
        if self.return_velocity:
            df[self.velocity_records] = u.velocities

        # TODO
        # S.S. to be implemented here !
        return df

        

    def get_topology_df(self) -> pd.DataFrame:
        df = {}
        if "name" in self.topology_records:
            df["name"] = self.atoms.names

        if self.return_alt:
            df["alt"] = self.atoms.altLocs
        
        if "resn" in self.topology_records:
            df["resn"] = self.atoms.resnames

        if self.return_chain:
            df["chain"] = self.atoms.chainIDs

        if "resi" in self.topology_records:
            df["resi"] = self.atoms.resids

        if self.return_occupancy:
            df["occupancy"] = self.atoms.occupancies

        if self.return_b:
            df["b"] = self.atoms.tempFactors

        if self.return_segi:
            df["segi"] = self.atoms.segids

        if self.return_e:
            df["e"] = self.atoms.elements

        if self.return_q:
            df["q"] = self.atoms.charges

        if self.return_m:
            df["m"] = self.atoms.masses


        df = pd.DataFrame(df)

        if "record_name" in self.topology_records:
            df["record_name"] = "ATOM"
            protein_resn = list(ONE_LETTER_CODE.keys())
            idx = df.query("not resn in @protein_resn").index
            df.loc[idx, "record_name"] = "HETATM"

            cols = list(df.columns)
            df = df[cols[-1:] + cols[:-1]]

        df.index += 1
        df.index.name = "atom_id"
        return df
        
    def reset(self, **kwargs):
        params = self.default_params.copy()
        for key, param in kwargs.items():
            if key in params:
                params[key] = param

        for key, param in params.items():
            setattr(self, key, param)

    def close(self):
        self.universe.trajectory.close()

    def load(self):
        self = iter(self)
        df = next(self)
        self.frame.close()

        return df
    
    def save(self, filename:str):
        raise NotImplementedError("`save` method not implemented for general Traj object. Please use one of the existing file format...")