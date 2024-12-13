import pandas as pd

try :
    from typing import Self

except ImportError:
    from typing_extensions import Self # python versions < 3.11

from io import TextIOWrapper

__all__ = ["Traj", "XYZ", "PDB", "GRO"]

class Traj:
    """
    A class to represent a molecular dynamics trajectory file handler.

    This class provides methods to open, read, and iterate over trajectory 
    files commonly used in molecular dynamics simulations.
    """
    def __init__(self, filename:str, mode = "r") -> None:
        """
        Initializes the Traj object and closes the file immediately.

        Args:
            filename (str): The name of the file containing trajectory data.
            mode (str, optional): The mode in which the file is opened. 
                Defaults to "r" (read mode).

        Returns:
            None
        """
        self.file = open(filename, mode)
        self.file.close()

    def __repr__(self) -> str:
        """
        Returns a string representation of the Traj object.

        The representation includes the class name, file name, and the file mode.

        Returns:
            str: A string representing the Traj object.
        """
        return f"md_manager.{self.__class__.__name__}: filename = '{self.file.name}', mode = '{self.file.mode}'"

    def __iter__(self) -> Self:
        """
        Makes the Traj object an iterable.

        Ensures that the file is open before iteration begins.

        Returns:
            Traj: The iterable object itself.
        """
        if self.file.closed:
            self.open()
        return self
    
    def __next__(self) -> pd.DataFrame:
        """
        Retrieves the next frame from the trajectory file.

        This method reads the next frame from the file and returns it as a 
        Pandas DataFrame. If the file is exhausted, raises a StopIteration 
        exception.

        Returns:
            pd.DataFrame: A DataFrame containing the data for the next frame.

        Raises:
            StopIteration: If there are no more frames in the trajectory file.
        """
        df =  self.next_frame(self.file)
        if len(df) > 0:
            return df
        
        # At this point the traj is over.
        raise StopIteration
    
    def open(self, mode="r"):
        """
        Opens the trajectory file.

        If the file is already open, this method will reopen it with the 
        specified mode.

        Args:
            mode (str, optional): The mode in which to open the file. 
                Defaults to "r" (read mode).

        Returns:
            None
        """
        self.file = open(self.file.name, mode)

    def close(self):
        """
        Closes the trajectory file.

        Ensures that the file resource is properly released.

        Returns:
            None
        """
        self.file.close()


######################################################################################################################################

class XYZ(Traj):
    """
    A specialized class for handling XYZ molecular dynamics trajectory files.

    Extends the Traj class to provide specific functionality for reading and 
    writing XYZ file formats, which store atomic coordinates and metadata.
    """
    columns = ["name", "x", "y", "z"]

    def __len__(self) -> int:
        """
        Counts the number of frames in an XYZ trajectory file.

        Returns:
            int: The number of frames in the trajectory.
        """
        if self.file.closed:
            self.open("r")
        frame_count = 0
        while True:
            try:
                line = self.file.readline()
                if not line:  # EOF
                    break
                Natom = int(line.strip())  # Read number of atoms
                for _ in range(Natom + 1):  # Skip atomic lines and comment line
                    _ = self.file.readline()
                frame_count += 1
            except ValueError:
                break
        self.file.seek(0)  # Reset file pointer to the beginning
        return frame_count

    @staticmethod
    def read_format(line:str) -> tuple[str, float, float, float]:
        """
        Parses a line from an XYZ file to extract atomic data.

        Args:
            line (str): A single line from the XYZ file containing atomic 
                information in the format: "name x y z".

        Returns:
            tuple[str, float, float, float]: A tuple containing the atom name 
            and its x, y, and z coordinates.
        """
        line = line.split()
        return (
            line[0],
            float(line[1]),
            float(line[1]),
            float(line[2])
        )
    
    @staticmethod
    def write_format(s:pd.Series) -> str:
        """
        Formats a Pandas Series representing an atom into an XYZ file line.

        Args:
            s (pd.Series): A Series containing the atom data with columns 
                "name", "x", "y", and "z".

        Returns:
            str: A formatted string to be written to an XYZ file.
        """
        return "".join([
            f"{s['name']:4s} ",
            f"{s['x']} ",
            f"{s['y']} ",
            f"{s['z']}\n"
        ])
    
    @classmethod
    def next_frame(cls, lines:TextIOWrapper|list[str]) -> pd.DataFrame:
        """
        Reads the next frame from an XYZ trajectory file.

        This method processes a set of lines from an XYZ file to extract atomic
        data into a Pandas DataFrame.

        Args:
            lines (TextIOWrapper | list[str]): The file or list of lines to 
                process. If a list is provided, it will be iterated over.

        Returns:
            pd.DataFrame: A DataFrame containing atomic data with columns 
            "name", "x", "y", and "z", and indexed by atom ID.

        Raises:
            StopIteration: If the input lines are exhausted before processing 
                a frame.
        """
        if type(lines) == list:
            lines = iter(lines)

        try :
            line = next(lines)
            Natom = int(line.strip())
            _ = next(lines)

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
        """
        Writes a single frame to the XYZ trajectory file.

        Args:
            df (pd.DataFrame): A DataFrame containing atomic data with columns 
                "name", "x", "y", and "z".
            model_id (int, optional): The model ID or timestep to associate 
                with the frame. Defaults to 1.

        Returns:
            None
        """
        if self.file.closed:
            self.open("w")

        self.file.write(f"{len(df)}\n")
        self.file.write(f"Generated by MD-manager. t = {model_id}\n")

        for _, atom in df.iterrows():
            self.file.write(XYZ.write_format(atom))




######################################################################################################################################

class PDB(Traj):
    """
    A specialized class for handling PDB (Protein Data Bank) trajectory files.

    This class provides methods for reading and writing atomic and molecular
    data in the PDB format, including support for standard fields such as 
    atom coordinates, chain IDs, and occupancy.
    """
    columns = ["record_name", "name", "alt", "resn", "chain", "resi", "insertion", "x", "y", "z", "occupancy", "b", "segi", "e", "q"]

    def __len__(self) -> int:
        """
        Counts the number of frames in a PDB trajectory file.

        Returns:
            int: The number of frames in the trajectory.
        """
        if self.file.closed:
            self.open("r")
        frame_count = 0
        for line in self.file:
            if line.startswith("MODEL"):
                frame_count += 1
        self.file.seek(0)  # Reset file pointer to the beginning
        return frame_count

    @staticmethod
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
        """
        Formats a Pandas Series representing an atom into a PDB file line.

        Args:
            id (int): The atom ID, typically a unique identifier within the 
                frame.
            s (pd.Series): A Series containing atomic data with columns 
                matching the `columns` attribute.

        Returns:
            str: A formatted string adhering to the PDB file format.
        """
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
        Reads atomic data from lines and returns the associated DataFrame.

        This method processes lines from a PDB file to extract atomic and 
        heteroatomic information until encountering an "END" record.

        Args:
            lines (TextIOWrapper | list[str]): The input source to process, 
                either as a file object or a list of strings.

        Returns:
            pd.DataFrame: A DataFrame containing atomic data with columns 
            matching `columns`. The index is labeled with atom IDs.
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
        Writes atomic data from a DataFrame into the PDB file format.

        Args:
            df (pd.DataFrame): A DataFrame containing atomic data with columns 
                matching `columns`.
            model_id (int, optional): The model ID to write into the PDB file, 
                typically corresponding to a frame or timestep. Defaults to 1.

        Returns:
            None
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
                new_line = PDB.write_format(atom_id, atom)
                self.file.write(new_line)
            self.file.write("TER\n")

        # Hetero atoms:
        for _, hetatm in HET.iterrows():
            atom_id += 1
            new_line = PDB.write_format(atom_id, hetatm)
            self.file.write(new_line)
        self.file.write("ENDMDL\n")

######################################################################################################################################

class GRO(Traj):
    """
    A specialized class for handling GRO (GROMACS) trajectory files.

    This class provides methods for reading and writing atomic data in the GRO 
    file format, which includes support for atomic coordinates, velocities, 
    and associated metadata.
    """
    columns = ["resi", "resn", "name", "x", "y", "z", "vx", "vy", "vz"]

    def __len__(self) -> int:
        """
        Counts the number of frames in a GRO trajectory file.

        Returns:
            int: The number of frames in the trajectory.
        """
        if self.file.closed:
            self.open("r")
        frame_count = 0
        while True:
            try:
                _ = self.file.readline()  # Skip header line
                if not self.file.readline():  # Check for empty atom count line
                    break
                Natom = int(self.file.readline().strip())  # Atom count
                for _ in range(Natom + 1):  # Skip atom lines and box line
                    _ = self.file.readline()
                frame_count += 1
            except ValueError:
                break
        self.file.seek(0)  # Reset file pointer to the beginning
        return frame_count

    @staticmethod
    def read_format(line: str) -> tuple[int, str, str, float, float, float, float, float, float]:
        """
        Parses a line from a GRO file to extract atomic data.

        Args:
            line (str): A single line from the GRO file, formatted according to 
                the GRO specification.

        Returns:
            tuple[int, str, str, float, float, float, float, float, float]: 
            A tuple containing:
                - Residue index (int)
                - Residue name (str)
                - Atom name (str)
                - X, Y, Z coordinates (float)
                - X, Y, Z velocities (float), defaulting to 0.0 if absent.
        """
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
        """
        Formats a Pandas Series representing an atom into a GRO file line.

        Args:
            id (int): The atom ID, typically a unique identifier within the frame.
            s (pd.Series): A Series containing atomic data with columns matching 
                the `columns` attribute.

        Returns:
            str: A formatted string adhering to the GRO file format.
        """
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
        Reads atomic data from lines and returns the associated DataFrame.

        This method processes lines from a GRO file to extract atomic 
        coordinates and velocities for a single frame.

        Args:
            lines (TextIOWrapper | list[str]): The input source to process, 
                either as a file object or a list of strings.

        Returns:
            pd.DataFrame: A DataFrame containing atomic data with columns 
            matching `columns`. The index is labeled with atom IDs.
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
        """
        Writes atomic data from a DataFrame into the GRO file format.

        This method writes a single frame to the GRO file, including atomic 
        coordinates and velocities.

        Args:
            df (pd.DataFrame): A DataFrame containing atomic data with columns 
                matching `columns`.
            model_id (int, optional): The model ID to write into the GRO file, 
                typically corresponding to a frame or timestep. Defaults to 1.

        Returns:
            None
        """
        if self.file.closed:
            self.open("w")

        self.file.write(f"Generated by MD-manager. t = {model_id}\n")
        self.file.write(f"{len(df)}\n")
        id = 0
        for _, atom in df.iterrows():
            id += 1
            self.file.write(GRO.write_format(id, atom))