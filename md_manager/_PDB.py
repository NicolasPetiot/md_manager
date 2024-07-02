from ._params import DF_COLUMNS, DF_TYPES, ATOMIC_MASSES

import pandas as pd
from io import TextIOWrapper

import sys
from urllib.request import urlopen

__all__ = [
    "PDB",
    "pdb2df",
    "df2pdb",
    "fetch_PDB",
    "COM",
    "shift_df",
    #rotate_df
]

class PDB:
    """class that contains methods for pdb file rreading and writing"""
    __write_format = {
        "record_name" : lambda x : f"{x:<6s}",
        "atom_id"     : lambda x : f"{x:5d}",
        "blank1"      : lambda _ : " ",
        "name"        : lambda x : f" {x:<3s}" if len(x) < 4 else f"{x:4s}",
        "alt"         : lambda x : f"{x:1s}",
        "resn"        : lambda x : f"{x:3s}",
        "blank2"      : lambda _ : " ",
        "chain"       : lambda x : f"{x:1s}",
        "resi"        : lambda x : f"{x:4d}",
        "insertion"   : lambda x : f"{x:1s}",
        "blank3"      : lambda _ : 3*" ",
        "x"           : lambda x : f"{x:8.3f}",
        "y"           : lambda x : f"{x:8.3f}",
        "z"           : lambda x : f"{x:8.3f}",
        "occupancy"   : lambda x : f"{x:6.2f}",
        "b"           : lambda x : f"{x:6.2f}",
        "blank4"      : lambda _ : 7*" ",
        "segi"        : lambda x : f"{x:<3s}",
        "e"           : lambda x : f"{x:2s}",
        "q"           : lambda x : f"{x:2s}"
    }
    
    def __init__(self, filename:str, mode = "r") -> None:
        self.file = open(filename, mode)
        self.file.close()

        self.__helix = []
        self.__sheet = []

    def __iter__(self):
        """called with iter(self)"""
        self.open()
        # initialise ss:
        helix, sheet = self.build_ss(line_iterator=self.file)
        self.__helix = helix
        self.__sheet = sheet
        return self
    
    def __next__(self):
        """called with next(self)"""
        if self.file.closed:
            iter(self)

        df = self.build_model(self.file, helix=self.__helix, sheet=self.__sheet)
        if len (df) > 0:
            return df
        
        else :
            self.close()
            raise StopIteration

    def open(self, mode = "r"):
        """set the file in open mode"""
        filename = self.file.name
        self.file = open(filename, mode)

    def close(self):
        self.file.close()

    ####################################################
    #                                                  #
    # Reading block: methods allowing PDB file reading #
    #                                                  #
    ####################################################

    @classmethod
    def build_model(cls, line_iterator:TextIOWrapper|list[str], helix = [], sheet = []) -> pd.DataFrame:
        """Creates DataFrame associated to the next frame in the pdb line iterator"""
        # Reading Data:
        atoms = []
        for line in line_iterator:
            record = line[:6].strip()

            # coord
            if record in {"ATOM", "HETATM"}:
                atoms.append(cls.read_coord(line, record))

        
            if record.startswith("END"):
                break
            
        # Building DataFrame:
        df = pd.DataFrame(
            atoms, 
            columns=DF_COLUMNS,
            index=pd.Index([1+i for i in range(len(atoms))], name="atom_id")
        ).astype(DF_TYPES)

        # Secondary structures:
        query_str = f"{DF_COLUMNS[4]} == @chain and {DF_COLUMNS[5]} >= @resiMin and {DF_COLUMNS[5]} <= @resiMax"
        if len(helix) > 0:
            s = pd.Series(False, index = df.index, name="helix")
            for chain, resiMin, resiMax in helix:
                # chain resiMin&Max being used by the query string
                idx = df.query(query_str).index
                s[idx] = True
            df["helix"] = s

        if len(sheet) > 0:
            s = pd.Series(False, index = df.index, name = "sheet")
            for chain, resiMin, resiMax in sheet:
                # chain resiMin&Max being used by the query string
                idx = df.query(query_str).index
                s[idx] = True
            df["sheet"] = s

        return df


    @staticmethod
    def read_coord(line:str, record:str) -> tuple:
        """
        Returns a tuple that contains the (record_name, name, alt, resn, chain, resi, insertion, x, y, z, occupancy, b, segi, elem, charge, mass) 
        informations about an atom line in the PDB file.

        By default,  the occupancy is set to 1.0.

        By default, the Bfactors are set to 0.0 AA².

        The mass is extracted from the `ATOMIC_MASSES` dictionary (see _params.py)
        
        args : 

        - `line:str` : Line of a pdb file that starts with "ATOM"/"HETATM"
        """
        #atom_id   = line[ 6:11]
        name      = line[12:16].strip()
        alt       = line[16:17].strip()
        resn      = line[17:20].strip()
        chain     = line[21:22].strip()
        resi      = int(line[22:26])
        insertion = line[26:27].strip()
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
        segi      = line[72:76].strip()
        elem      = line[76:78].strip()
        charge    = line[78:80].strip()
        try :
            mass      = ATOMIC_MASSES[elem]
        except KeyError:
            raise KeyError(f"Unknown element symbol '{elem}'. Please update the ATOMIC_MASSES dictionary in '_parameters.py'.")
        return record, name, alt, resn, chain, resi, insertion, x, y, z, occupancy, b, segi, elem, charge, mass
    
    @classmethod
    def build_ss(cls, line_iterator:TextIOWrapper|list[str]) -> tuple[list[str], list[str]]:
        """Creates lists describing the position of helices and sheets in the line iterator."""
        helix = []
        sheet = []
        for i, line in enumerate(line_iterator):
            record = line[:6].strip()

            if record == "HELIX":
                helix.append(cls.read_helix(line))

            if record == "SHEET":
                sheet.append(cls.read_sheet(line))

            if record == "MODEL":
                break
            
            if record == "ATOM" and type(line_iterator) == TextIOWrapper:
                line_iterator.seek(i-1) # back to the previous line!
                break

        return helix, sheet
    
    @staticmethod
    def read_helix(line:str) -> tuple[str, int, int]:
        """Returns relevant informations encoded in HELIX pdb line"""

        chainMin = line[19]
        resiMin = int(line[21:25])

        chainMax = line[31]
        resiMax = int(line[33:37])

        if chainMin != chainMax:
            raise ValueError("Helix in differrent chains")
        
        return chainMin, resiMin, resiMax
    
    @staticmethod
    def read_sheet(line:str) -> tuple[str, int, int]:
        """return relevant informations encoded in SHEET pdb line"""
        chainMin = line[21]
        resiMin = int(line[22:26])

        chainMax = line[32]
        resiMax = int(line[33:37])

        if chainMin != chainMax:
            raise ValueError("Helix in differrent chains")
        
        return chainMin, resiMin, resiMax

    ####################################################
    #                                                  #
    # Writing block: methods allowing PDB file writing #
    #                                                  #
    ####################################################

    def write_pdb(self, df:pd.DataFrame=None, dfs:list[pd.DataFrame]=None):
        """Creates a structure/trajectory file in pdb format."""
        self.open("w")
        is_traj = False
        if type(dfs) == list:
            is_traj = True
            df = dfs[0]

        for line in self.build_ss_lines(df):
            self.file.write(line)

        for line in self.build_model_lines(df):
            self.file.write(line)

        if is_traj:
            for i, df in enumerate(dfs[1:]):
                for line in self.build_model_lines(df, model_id=i+2):
                    self.file.write(line)
        self.close()

    @classmethod
    def build_ss_lines(cls, df:pd.DataFrame)->list[str]:
        """Returs a list of line to be written in a pdb file to describe secodary structures"""
        lines = []
        if "helix" in df.columns:
            query_str = f"{DF_COLUMNS[1]} == 'CA' and helix"
            col_slice = DF_COLUMNS[3:6]

            ss = df.query(query_str)[col_slice]
            resMin = ss[ss.resi.diff()   !=  1].values # starting points of helices
            resMax = ss[ss.resi.diff(-1) != -1].values # end point of helices

            for i, ((resnMin, chainMin, resiMin), (resnMax, chainMax, resiMax)) in enumerate(zip(resMin, resMax)):
                line = f"HELIX  {i+1:3d} {cls.ss_identifier_string(i)}"

                # initial residue:
                line += f" {resnMin:3s} {chainMin:s} {resiMin:4d}"

                # last residue:
                line += f"  {resnMax:3s} {chainMax:s} {resiMax:4d}"
                
                line += "  1" # kind of helix not in memory (1 by default)

                length = resiMax - resiMin
                line += 31*" "
                line += f"{length+1:5d}"
                lines.append(line + "\n")

        if "sheet" in df.columns:
            query_str = f"{DF_COLUMNS[1]} == 'CA' and sheet"
            col_slice = DF_COLUMNS[3:6]

            ss = df.query(query_str)[col_slice]
            resMin = ss[ss.resi.diff()   !=  1].values # starting points of sheets
            resMax = ss[ss.resi.diff(-1) != -1].values # end point of sheets

            for i,((resnMin, chainMin, resiMin), (resnMax, chainMax, resiMax)) in enumerate(zip(resMin, resMax)):
                print(resnMin, chainMin, resiMin)
                line = f"SHEET  {i+1:3d} {cls.ss_identifier_string(i)} 1"

                # initial residue:
                line += f" {resnMin:3s} {chainMin:s}{resiMin:4d}"

                # last residue:
                line += f"  {resnMax:3s} {chainMax:s}{resiMax:4d}"

                line += "  0" # default sheet type
                lines.append(line + "\n")

        return lines
    
    @classmethod
    def build_model_lines(cls, df:pd.DataFrame, model_id = 1) -> list[str]:
        """Returns a list of lines to be written in a file to describe atoms coordiates."""
        df.index = pd.Index([1+i for i in range(len(df))], name = "atom_id") # make sure the atom idx are correct
        APO = df.query("record_name == 'ATOM'")
        HET = df.query("record_name == 'HETATM'")

        lines = [f"MODEL{model_id:8d}\n"]
        for _, chain in APO.groupby("chain"):
            for _, atom in chain.iterrows():
                lines.append(cls.__build_coord_line(atom))
            lines.append("TER\n")

        for _, atom in HET.iterrows():
            lines.append(cls.__build_coord_line(atom))
        lines.append("ENDMDL\n")
        return lines
    
    @classmethod
    def __build_coord_line(cls, atom:pd.Series) -> str:
        """Returns a string describing an atom in pdb format."""
        line = ""
        for col, format in cls.__write_format.items():
            if col.startswith("blank"):
                line += format(col)

            elif col.startswith("atom_id"):
                line += format(atom.name)

            else:
                line += format(atom[col])

        return line + "\n"

    @staticmethod
    def ss_identifier_string(id:int)->str:
        """Converts a number into string identifier for pdb file"""
        n = id % 9
        remind = (id - n)//9

        b = remind % 26
        remind = (remind - b)//26

        a = remind % 26

        return f"{chr(65+a)}{chr(65+b)}{n+1}"
    
####################################################
#                                                  #
# Functions that make the use of MD-manager easier #
#                                                  #
####################################################

def pdb2df(filename:str, atom_only = False) -> pd.DataFrame:
    """Read all the atoms in the first MODEL of a pdb file"""
    pdb = PDB(filename)
    df = next(pdb)

    if atom_only:
        query_str = f"{DF_COLUMNS[0]} == 'ATOM'"
        df = df.query(query_str)

    return df

def df2pdb(filename:str, data:pd.DataFrame|list[pd.DataFrame]):
    """Generates a pdb structure/trajctory according to the type of the input data"""
    new = PDB(filename, "w")
    if type(data) == list:
        new.write_pdb(dfs=data)
    
    else:
        new.write_pdb(df=data)


def fetch_PDB(pdb_code:str, atom_only = False) -> pd.DataFrame:
    """
    Returns a pandas.DataFrame associated to the structure loaded from the 'rcsb.org' website.
    """
    url = f"https://files.rcsb.org/download/{pdb_code.lower()}.pdb"
    response = urlopen(url)

    txt = response.read()
    lines = (txt.decode("utf-8") if sys.version_info[0] >= 3 else txt.decode("ascii")).splitlines()

    helix, sheet = PDB.build_ss(lines)
    df = PDB.build_model(lines, helix, sheet)

    if atom_only:
        query_str = f"{DF_COLUMNS[0]} == 'ATOM'"
        df = df.query(query_str)

    return df

def COM(df:pd.DataFrame) -> pd.Series:
    """Returns the center of mass of an input DataFrame"""
    xyz = ["x", "y", "z"]
    com = pd.Series(0.0, index=xyz+["m"])

    for _, atom in df.iterrows():
        pos = atom[xyz].astype(float)
        m = atom.m

        com[xyz] += pos*m
        com.m += m

    com[xyz] /= com.m
    return com

def shift_df(df:pd.DataFrame, vec:pd.Series) -> pd.DataFrame:
    """"""
    xyz = ["x", "y", "z"]
    df[xyz] = df.apply(
        lambda s: s[xyz] + vec[xyz]
    )
    return df

#TODO def rotate_df(df:pd.DataFrame, (...)) -> pd.DataFrame: