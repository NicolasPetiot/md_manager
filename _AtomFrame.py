from ._parameters import *
import pandas as pd
import numpy as np

class PDBfile :
    """
    Iterator object that allow to loop over the frames in a md trajectory in the PDB format.

    PDBfile(filename).read2df() returns a frame that contains all the atoms of the file.
    """
    def __init__(self, filename) -> None:
        self.file = open(filename, "r")
        self.file.close()

    def __len__(self) -> int :
        """
        Returns the number of frames in trajectory.
        """
        len = 0 
        close = False
        if self.file.closed :
            self.file = open(self.file.name, "r")
            close = True

        for line in self.file :
            if line[:5] == "MODEL":
                len += 1

        if close :
            self.file.close()

        return len

    ### Allow to use 'with statments' ###
    def __enter__(self):
        """
        Called using 'with MdTraj(filename) as self'.
        """
        self.file = open(self.file.name, "r")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Called at the end of a 'with MdTraj(filename) as self' statment.
        """
        self.file.close()

    
    #####################################

    ### Allow to loop over frames ###
    def __iter__(self):
        """
        Is called at the beggining of a 'for model in pdb' statment.

        Caution : pdb.file must be oppened.
        """
        if self.file.closed:
            self.file = open(self.file.name, "r")
        return self
    
    def __next__(self):
        """
        Is called using 'for model in pdb'.
        Returns a pandas.DataFrame that corresponds to the next MODEL in PDB file.

        Caution : pdb.file must be oppened.

        Caution : the text readed must contains 'ENDMDL' at the end of a model. If not, please use 'read2df' instead.
        """
        atoms = []
        for line in self.file :
            if line[:6] == "ENDMDL":
                return self._build_df_from_atom_list(atoms)
            
            # append atoms with potential new informations
            if line[:6] in {"ATOM  ", "HETATM"}:
                atoms.append(self._scan_pdb_line(line))
            
        raise StopIteration
    
    def open(self):
        """
        Opens the pdb.file wrapper.
        """
        self.file = open(self.file.name, "r")

    def close(self):
        """
        Closes the pdb.file wrapper.
        """
        self.file.close()
    #################################

    ### Manipulate pdb lines ###
    @staticmethod
    def _build_df_from_atom_list(atoms:list):
        """
        Returns a pandas.DataFrame that contains the atoms in the input list.
        """
        columns=["record_name", "name", "alt", "resn", "chain", "resi", "insertion", "x", "y", "z", "occupancy", "b", "segi", "e", "q", "m"]
        return pd.DataFrame(atoms, columns=columns)

    @staticmethod
    def _scan_pdb_line(line:str):
        """
        Returns a tuple that contains the (record_name, name, alt, resn, chain, resi, insertion, x, y, z, occupancy, b, segi, elem, charge, mass) 
        informations about an atom line in the PDB file.
        """
        if line.startswith("ATOM"):
            record_name = "ATOM  "
        elif line.startswith("HETATM"):
            record_name = "HETATM"
        else :
            raise ValueError("Input 'line' is not associated to an atom in a pdb file.")
        
        #atom_id   = line[ 6:11]
        name      = line[12:16]
        alt       = line[16:17]
        resn      = line[17:20]
        chain     = line[21:22]
        resi      = int(line[22:26])
        insertion = line[26:27]
        x         = float(line[30:38])
        y         = float(line[38:46])
        z         = float(line[46:54])
        try :
            occupancy = float(line[54:60])
        except ValueError :
            occupancy = 1.0
        try :
            b         = float(line[60:66])
        except ValueError :
            b         = 0.0
        segi      = line[72:76]
        elem      = line[76:78]
        charge    = line[78:79]
        try :
            mass      = ATOMIC_MASSES[elem]
        except KeyError:
            raise KeyError(f"Unknown element symbol '{elem}'. Please update the ATOMIC_MASSES dictionary in '_parameters.py'.")
        return record_name, name, alt, resn, chain, resi, insertion, x, y, z, occupancy, b, segi, elem, charge, mass
    
    def read2df(self):
        """
        Returns a pandas.DataFrame that contains all the atoms fount in the PDB file.
        """
        atoms = []
        with open(self.file.name, "r") as file :
            for line in file:
                if line[:6] in {"ATOM  ", "HETATM"}:
                    atoms.append(self._scan_pdb_line(line))
        return self._build_df_from_atom_list(atoms) 
    
    @classmethod
    def df2pdb(cls, df:pd.DataFrame, filename: str, title = "", remark = "") -> None:
        """
        Generates a pdb file from an inputed DataFrame.
        """
        with open(filename, "w") as file :
            if title != "":
                file.write("TITLE     " + title + "\n")
            if remark != "":
                for line in remark.splitlines():
                    file.write("REMARK    " + line + "\n")
            
            file.write("MODEL        1\n")
            id = 0
            APO = df[df.record_name == "ATOM  "]
            HET = df[df.record_name == "HETATM"]
            for _, chain in APO.groupby(["chain"]):
                for _, atom in chain.iterrows():
                    id += 1
                    line = cls._generate_pdb_line(atom, id)
                    file.write(line)
                file.write("TER\n")
            for atom in HET.iterrows():
                id += 1
                line = cls._generate_pdb_line(atom, id)
                file.write(line)
            file.write("TER\n")
            file.write("ENDMDL\n")

    @staticmethod
    def _generate_pdb_line(atom:pd.Series, id:int) -> str:
        """
        Generates a pdb line from an input series that contains atom's information.
        """
        # extract strings :
        record_name = atom["record_name"]
        name        = atom["name"]
        alt         = atom["alt"]
        resn        = atom["resn"]
        chain       = atom["chain"]
        insertion   = atom["insertion"]
        segi        = atom["segi"]
        elem        = atom["e"]
        charge      = atom["q"]

        # extract float/int
        id    = f"{id:5d}"
        resi  = f"{atom.resi:4d}"
        x     = f"{atom.x:8.3f}"
        y     = f"{atom.y:8.3f}"
        z     = f"{atom.z:8.3f}"
        occup = f"{atom.occupancy:6.2f}"
        b     = f"{atom.b:6.2f}"

        line = "%s%s %s%s%s %s%s%s   %s%s%s%s%s      %s%s%s\n"%(
            record_name, id, name, alt, resn, chain, resi, insertion, x, y, z, occup, b, segi, elem, charge
        )
        return line


    ############################
    

def atom_position(df:pd.DataFrame) -> np.ndarray:
    """
    Returns a numpy.ndarray that contains the [x, y, z] coordinates of the atoms in df.
    """
    return df[["x", "y", "z"]].to_numpy(dtype=float)

def atom_mass(df:pd.DataFrame) -> np.ndarray:
    """
    Returns a numpy.ndarray that contains the masses of the atoms in df.
    """
    return df["m"].to_numpy(dtype=float)

def atom_bfactors(df:pd.DataFrame) -> np.ndarray:
    """
    Returns a numpy.ndarray that contains the thermal bfactors of the atoms in df.
    """
    return df["b"].to_numpy(dtype=float)

    