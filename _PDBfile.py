from .utils._PDBfile_utils import _build_df_from_atom_list, _scan_pdb_line, _generate_pdb_line_list

__all__ = ["pdb"]

class pdb:
    """
    
    """
    def __init__(self, filename:str, write = False) -> None:
        try :
            self.file = open(filename, "r")
            self.file.close()
        
        except FileNotFoundError:
            if write:
                self.file = open(filename, "w")
                self.file.close()
            else :
                raise FileNotFoundError("File doesn't exists. Please set `write = True` to create a new file.")
        

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
    
    def __enter__(self):
        """
        Called using `with md.pdb as pdb:` statment...
        """
        self.file = open(self.file.name, "r")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Called at the end of a `with md.pdb as pdb:` statment...
        """
        self.file.close()

    def __next__(self):
        """
        Called using a `for model in pdb` statment.

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

    def read2df(self):
        """
        
        """
        if self.file.closed:
            raise ValueError("Please use `open(self)` before reading file...")
        
        atoms = []
        for line in self.file:
            if line[:6] in {"ATOM  ", "HETATM"}:
                atoms.append(_scan_pdb_line(line))
        return _build_df_from_atom_list(atoms)
    

    def write_pdb(self, df=None, dfs=None, title=None):
        """
        
        """
        NoneType = type(None)
        header = "REMARK Generated by md_manager\n"
        lines = []
        if type(df) != NoneType:
            lines += _generate_pdb_line_list(df)

        elif type(dfs) != NoneType:
            for i, df in enumerate(dfs) :
                lines += ["MODEL %6d\n"%(i+1)]
                lines += _generate_pdb_line_list(df)
                lines += ["ENDMDL\n"]

        else :
            raise ValueError("No input DataFrame(s)")

        with open(self.file.name, "w") as file :
            if title != None:
                file.write(title)
            file.write(header)
            for line in lines:
                file.write(line)
