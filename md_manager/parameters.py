__all__ = [
    "TEMPERATURE",
    "BOLTZMANN",
    "AVOGADRO",
    "ONE_LETTER_CODE",
    "ATOM_NAME_CHI"
]

TEMPERATURE = 300.0          # K
BOLTZMANN   = 1.987204259e-3 # kcal/mol/K (see : 'https://en.wikipedia.org/wiki/Boltzmann_constant')
AVOGADRO    = 6.02214076e23  # 1/mol      (see : 'https://en.wikipedia.org/wiki/Avogadro_constant')

ONE_LETTER_CODE ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q',
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',
'GLY':'G', 'PRO':'P', 'CYS':'C'}

ATOM_NAME_CHI = ["N", "CA", "CB", "CG", "SG", "CG1", "OG", "OG1", "CD", "SD", "CD1", "OD1", "ND1", "CE", "NE", "OE1", "CZ", "NZ", "NH1"]

# obselete
ATOM_NAME_SELECTION_CHI = {
    "ALA" : ["CA"],
    "ARG" : ["N", "CA", "CB", "CG", "CD", "NE", "CZ", "NH1"],
    "ASN" : ["N", "CA", "CB", "CG", "OD1"],
    "ASP" : ["N", "CA", "CB", "CG", "OD1"],
    "CYS" : ["N", "CA", "CB", "SG"],
    "GLN" : ["N", "CA", "CB", "CG", "CD", "OE1"],
    "GLU" : ["N", "CA", "CB", "CG", "CD", "OE1"],
    "GLY" : ["CA"],
    "HIS" : ["N", "CA", "CB", "CG", "ND1"],
    "ILE" : ["N", "CA", "CB", "CG1", "CD"],
    "LEU" : ["N", "CA", "CB", "CG", "CD1"],
    "LYS" : ["N", "CA", "CB", "CG", "CD", "CE", "NZ"],
    "MET" : ["N", "CA", "CB", "CG", "SD", "CE"],
    "PHE" : ["N", "CA", "CB", "CG", "CD1"],
    "PRO" : ["N", "CA", "CB", "CG", "CD"],
    "SER" : ["N", "CA", "CB", "OG"],
    "THR" : ["N", "CA", "CB", "OG1"],
    "TRP" : ["N", "CA", "CB", "CG", "CD1"],
    "TYR" : ["N", "CA", "CB", "CG", "CD1"],
    "VAL" : ["N", "CA", "CB", "CG1"]
} # atoms used for computation of dihedral angles (see 'http://www.mlb.co.jp/linux/science/garlic/doc/commands/dihedrals.html')

import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")