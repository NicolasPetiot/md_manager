import pandas as pd

from typing import Self
from io import TextIOWrapper

__all__ = ["XYZ", "PDB", "GRO"]

class Traj:
    def __init__(self, filename:str, mode = "r") -> None:
        self.file = open(filename, mode)
        self.file.close()

    def __repr__(self) -> str:
        return f"md_manager.{self.__class__.__name__}: filename = '{self.file.name}', mode = '{self.file.mode}'"

    def __iter__(self) -> Self:
        """
        Must ensure that Self is an open file.
        """
        if self.file.closed:
            self.open()
        return self
    
    def __next__(self) -> pd.DataFrame:
        df =  self.next_frame(self.file)
        if len(df) > 0:
            return df
        
        # At this point the traj is over.
        raise StopIteration
    
    def open(self, mode="r"):
        self.file = open(self.file.name, mode)

    def close(self):
        self.file.close()


######################################################################################################################################

class XYZ(Traj):
    columns = ["name", "x", "y", "z"]

    @staticmethod
    def read_format(line:str) -> tuple[str, float, float, float]:
        """Read atom line in XYZ file"""
        line = line.split()
        return (
            line[0],
            float(line[1]),
            float(line[1]),
            float(line[2])
        )
    
    @staticmethod
    def write_format(_:int, s:pd.Series) -> str:
        """Create a line to be written in an XYZ file"""
        return "".join([
            f"{s["name"]:5s} ",
            f"{s["x"]:12.6f} ",
            f"{s["y"]:12.6f} ",
            f"{s["z"]:12.6f}\n"
        ])
    
    @classmethod
    def next_frame(cls, lines:TextIOWrapper|list[str]) -> pd.DataFrame:
        #TODO
        raise NotImplemented("Must implement XYZ reading")

    def write_frame(self, df:pd.DataFrame) -> None:
        #TODO
        raise NotImplemented("Must implement XYZ writing")

######################################################################################################################################

class PDB(Traj):
    columns = ["record_name", "name", "alt", "resn", "chain", "resi", "insertion", "x", "y", "z", "occupancy", "b", "segi", "e", "q"]

    @staticmethod
    def read_format(line:str) -> tuple[str, str, str, str, int, str, float, float, float, float, float, float, str, str, str]:
        """Read atom line in PDB file"""
        return (
            #line[:6].strip(),
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
    
    @staticmethod
    def write_format(id:int, s:pd.Series) -> str:
        """Creates a line in PDB format"""
        return "".join([
            f"{s['record_name']:<6s}",   # record_name (ATOM/HETATM)
            f"{id:>5d} ",                # atom_id
            f" {s['name']:<3s}" if len(s['name']) < 4 else f"{s['name']:4s}",  # atom_name
            f"{s['alt']:1s}",            # alt_loc
            f"{s['resn']:<3s} ",         # resn (residue name)
            f"{s['chain']:1s}",          # chain_id
            f"{s['resi']:>4d}",          # resi (residue sequence number)
            f"{s['insertion']:1s}   ",   # insertion code (often blank)
            f"{s['x']:>8.3f}",           # x coordinate
            f"{s['y']:>8.3f}",           # y coordinate
            f"{s['z']:>8.3f}",           # z coordinate
            f"{s['occupancy']:>6.2f}",   # occupancy
            f"{s['b']:>6.2f}          ", # b-factor (temp factor)
            f"{s['e']:>2s}",             # element symbol
            f"{s['q']:>2s}\n"            # charge
        ])

    @classmethod
    def next_frame(cls, lines:TextIOWrapper|list[str]) -> pd.DataFrame:
        """
        Read the atoms in lines and returns the associated DataFrame.
        """

        atoms = []
        for line in lines:
            record = line[:6].strip()
            if record.startswith("END"):
                break

            elif record in {'ATOM', 'HETATM'}:
                name, alt, resn, chain, resi, insertion, x, y, z, occupancy, b, segi, e, q = cls.read_format(line)
                atoms.append([record, name, alt, resn, chain, resi, insertion, x, y, z, occupancy, b, segi, e, q])

        df = pd.DataFrame(atoms, columns=cls.columns)
        # TODO: add S.S.

        df.index = pd.Index([i for i in range(1, len(df)+1)], name = "atom_id")
        return df
    
    def write_frame(self, df:pd.DataFrame, model_id = 1):
        """
        
        """
        if self.file.closed:
            self.open(mode="w")
            self.file.write("REMARK generated by MD-manager\n")

        self.file.write("MODEL %8d\n"%model_id)

        APO = df.query("record_name == 'ATOM'")
        HET = df.query("record_name == 'HETATM'")

        atom_id = 0
        # Atoms -> must split chains:
        for _, group in APO.groupby("chain"):
            for _, atom in group.iterrows():
                atom_id += 1
                new_line = PDB.write_format(id, atom)
                self.file.write(new_line)
            self.file.write("TER\n")

        # Hetero atoms:
        for _, hetatm in HET.iterrows():
            atom_id += 1
            new_line = PDB.write_format(id, hetatm)
            self.file.write(new_line)
        self.file.write("ENDMDL\n")

######################################################################################################################################

class GRO(Traj):
    columns = ["resi", "resn", "name", "x", "y", "z", "vx", "vy", "vz"]

    @staticmethod
    def read_format(line: str) -> tuple[int, str, str, float, float, float, float, float, float]:
        """Read atom line in GRO file"""
        return (
            int(line[0:5].strip()),     # residue index
            line[5:10].strip(),         # residue number
            line[10:15].strip(),        # atom name
            #int(line[15:20].strip()),   # atom id
            float(line[20:28].strip()), # x coordinate
            float(line[28:36].strip()), # y coordinate
            float(line[36:44].strip()), # z coordinate
            float(line[44:52].strip()) if line[44:52].strip() != "" else 0.0, # x velocity
            float(line[52:60].strip()) if line[52:60].strip() != "" else 0.0, # y velocity
            float(line[60:68].strip()) if line[60:68].strip() != "" else 0.0  # z velocity
        )

    @staticmethod
    def write_format(id:int, s:pd.Series) -> str:
        """Creates line for GRO file"""
        return "".join([
            f"{s['resi']:>5d}", # Residue number
            f"{s['resn']:<5s}", # Residue name
            f"{s['name']:>5s}", # Atom name
            f"{id:>5d}",        # Atom number
            f"{s['x']:8.3f}",   # X coordinate
            f"{s['y']:8.3f}",   # Y coordinate
            f"{s['z']:8.3f}",   # Z coordinate
            f"{s['vx']:8.4f}" if 'vx' in s else "",  # X velocity, if present
            f"{s['vy']:8.4f}" if 'vy' in s else "",  # Y velocity, if present
            f"{s['vz']:8.4f}" if 'vz' in s else "",  # Z velocity, if present
            "\n"
        ])
    
    @classmethod
    def next_frame(cls, lines:TextIOWrapper|list[str]) -> pd.DataFrame:
        """
        Read the atoms in lines and returns the associated DataFrame.
        """
        if type(lines) == list:
            lines = iter(lines)

        try :
            _ = next(lines) # header
            line = next(lines)
            Natom = int(line.strip())

        except StopIteration:
            Natom = 0

        atoms = []
        for _ in range(Natom):
            line = next(lines)
            atoms.append(cls.read_format(line))

        df = pd.DataFrame(atoms, columns=cls.columns)
        df.index = pd.Index([i for i in range(1, len(df)+1)], name = "atom_id")
        return df
    
    def write_frame(self, df:pd.DataFrame, model_id = 1) -> None:
        """"""
        if self.file.closed:
            self.open("w")

        self.file.write(f"Generated by MD-manager. t = {model_id}\n")
        self.file.write(f"{len(df)}\n")
        id = 0
        for _, atom in df.iterrows():
            id += 1
            self.file.write(GRO.write_format(id, atom))