from ._parameters import *
from ._utils import *

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
                return _build_df_from_atom_list(atoms)
            
            # append atoms with potential new informations
            if line[:6] in {"ATOM  ", "HETATM"}:
                atoms.append(_scan_pdb_line(line))
            
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
    def read2df(self):
        """
        Returns a pandas.DataFrame that contains all the atoms fount in the PDB file.
        """
        atoms = []
        with open(self.file.name, "r") as file :
            for line in file:
                if line[:6] in {"ATOM  ", "HETATM"}:
                    atoms.append(_scan_pdb_line(line))
        return _build_df_from_atom_list(atoms) 
    ############################
    
