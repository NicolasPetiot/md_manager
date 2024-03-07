from .utils._params import DF_COLUMNS, DF_TYPES, ATOMIC_MASSES, STRING_FORMAT

from pandas import DataFrame, Series

class PDB:
    """
    
    """
    def __init__(self, filename:str, write=False) -> None:
        """
        Creates an instance of the `PDB` class.

        args :

        - `filename:str` : File name or relative path of the pdb file to read/write.

        - `write:bool = False` : Allow the creation of a the file is file not found.
        """
        # Open filename
        try :
            self.__file = open(filename, "r")
            self.__file.close()

        except FileNotFoundError as e:
            if write:
                self.__file = open(filename, "w")
                self.__file.close()

            else :
                raise e
        
        self.__is_open = False
        self.__returned_once = False

    def __iter__(self):
        """
        Initialises the iteration over models.
        """
        if not self.__is_open:
            self.open()
        self.__returned_once = False

        return self
    
    def __next__(self):
        """
        Iterates over the lines in the file and returns the DataFrame associated for each 'ENDMDL' in the file.
        """
        atoms = []
        
        for line in self.__file:
            if line.startswith("ENDMDL"):
                self.__returned_once = True
                return self.build_df_from_atoms(atoms)
            
            if line[:6] in {"ATOM  ", "HETATM"}:
                atoms.append(self.scan_pdb_line(line))

        else :
            if not self.__returned_once:
                self.__returned_once = True
                return self.build_df_from_atoms(atoms)
            self.close()
            raise StopIteration
            
    
    def open(self, mode = "r"):
        """
        Set the I/O wrapper associated to the instance of PDB in open mode.

        args : `mode:str = "r"` I/O interation mode (read by default)
        """
        self.__file = open(self.__file.name, mode)
        self.__is_open = True
    
    def close(self):
        """
        Set the I/O wrapper associated to the instance of PDB in close mode.
        """
        self.__file.close()
        self.__is_open = False
    
    def write(self, model = None, model_list = None):
        """
        Read the input DataFrame(s) and creates the associated pdb file.

        args : 

        -  `model:DataFrame` : Single frame to be written in the output file.

        - `model_list:list[DataFrame]` : List of frames to be written in the output file.
        """
        lines = []
        model_id = 1
        if model is not None:
            lines += self.generate_atom_lines(model, model_id)

        elif model_list is not None :
            for model in model_list:
                lines += self.generate_atom_lines(model, model_id)
                model_id += 1

        self.open("w")
        self.__file.writelines(lines)
        self.close()


    @staticmethod
    def build_df_from_atoms(atoms) -> DataFrame:
        """
        Creates the DataFrame associated to input `atoms` with expected column names and types.

        args :

        - `atoms:list[tuple]` : List created from `scan_pdb_line` method.
        """
        return DataFrame(atoms, columns=DF_COLUMNS).astype(DF_TYPES)
    
    @staticmethod
    def scan_pdb_line(line:str) -> tuple:
        """
        Returns a tuple that contains the (record_name, name, alt, resn, chain, resi, insertion, x, y, z, occupancy, b, segi, elem, charge, mass) 
        informations about an atom line in the PDB file.

        By default,  the occupancy is set to 1.0.

        By default, the Bfactors are set to 0.0 AA².

        The mass is extracted from the `ATOMIC_MASSES` dictionary (see _params.py)
        
        args : 

        - `line:str` : Line of a pdb file that starts with "ATOM"/"HETATM"
        """
        if line.startswith("ATOM"):
            record_name = "ATOM"
        elif line.startswith("HETATM"):
            record_name = "HETATM"
        else :
            raise ValueError("Input 'line' is not associated to an atom in a pdb file.")
        
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
        return record_name, name, alt, resn, chain, resi, insertion, x, y, z, occupancy, b, segi, elem, charge, mass
    
    @staticmethod
    def generate_atom_line(atom:Series, atom_id:int) -> str:
        """
        Generates a pdb line from an input series that contains atom's information.
        
        args : 

        - `atom:Series` : Series associated to a line of a DataFrame.

        - `atom_id:int` : Atom index in the pdb file.
        """
        line = ""
        col = DF_COLUMNS[0]
        line += STRING_FORMAT[col](atom[col])
        line += f"{atom_id:5d}"

        # blank 1:
        line += " "

        for col in DF_COLUMNS[1:4]:
            data = atom[col]
            line += STRING_FORMAT[col](data)

        # blank 2:
        line += " "
        
        for col in DF_COLUMNS[4:7]:
            data = atom[col]
            line += STRING_FORMAT[col](data)

        # blank 3:
        line += 3*" "

        for col in DF_COLUMNS[7:13]:
            data = atom[col]
            line += STRING_FORMAT[col](data)

        # blank 4:
        line += 7*" "

        for col in DF_COLUMNS[13:-1]:
            data = atom[col]
            line += STRING_FORMAT[col](data)

        return line + "\n"
    
    @classmethod
    def generate_atom_lines(cls, df:DataFrame, model_id:int = 1) -> list[str]:
        """
        Generates a list of lines to be written in a pdb file from input DataFrame and model index.

        args :

        - `df:DataFrame` : DataFrame to convert in lines in pdb format.

        - `model_id:int = 1` : Number of the model associated to the current frame.
        """
        lines = [f"MODEL{model_id:8d}\n"]
        # split atoms and hetero atoms :
        query_string = f"{DF_COLUMNS[0]} == 'ATOM'"
        APO = df.query(query_string)
        query_string = f"{DF_COLUMNS[0]} == 'HETATM'"
        HET = df.query(query_string)

        id = 0
        for _, chain in APO.groupby(DF_COLUMNS[4]):
            for _, atom in chain.iterrows():            
                id += 1
                lines.append(cls.__generate_atom_line(atom, atom_id=id))
            lines.append("TER\n")
        for _, atom in HET.iterrows():
            id += 1
            lines.append(cls.__generate_atom_line(atom, atom_id=id))
        lines.append("ENDMDL\n")
                
        return lines