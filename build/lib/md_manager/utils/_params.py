__all__ = [
    "TEMPERATURE",
    "BOLTZMANN",
    "AVOGADRO",
    "ATOMIC_MASSES",
    #"ATOM_NAME_SELECTION_CHI",
    "DF_COLUMNS",
    "DF_TYPES",
    "ONE_LETTER_CODE",
    #"THREE_LETTERS_CODE",
    "STRING_FORMAT",
]

TEMPERATURE = 300.0          # K
BOLTZMANN   = 1.987204259e-3 # kcal/mol/K (see : 'https://en.wikipedia.org/wiki/Boltzmann_constant')
AVOGADRO    = 6.02214076e23  # 1/mol      (see : 'https://en.wikipedia.org/wiki/Avogadro_constant')

ATOMIC_MASSES = {
    "H" :  1.0079,
    "C" : 12.0107,
    "N" : 14.0067,
    "O" : 15.9994,
    "S" : 32.0650,

    "NA" : 22.9897, # Sodium
    "MG" : 24.3050, # Magnesium
    "CL" : 35.4530, # Chlorine
    "K" : 39.0983, # Potassium
    "AR" : 39.9480, # Argon
    "CA" : 40.0780, # Calcium
    "CU" : 63.5460, # Copper

    "" : 1.0, 
} # g/mol (see 'https://www.lenntech.com/periodic/mass/atomic-mass.htm')

ONE_LETTER_CODE ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q',
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',
'GLY':'G', 'PRO':'P', 'CYS':'C'}

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
    "ILE" : ["N", "CA", "CB", "CG1", "CD1"],
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

DF_COLUMNS = ["record_name", "name", "alt", "resn", "chain", "resi", "insertion", "x", "y", "z", "occupancy", "b", "segi", "e", "q", "m"]
type_list   = [str, str, str, str, str, int, str, float, float, float, float, float, str, str, str, float]
DF_TYPES = {col:_type for col, _type in zip(DF_COLUMNS, type_list)}

STRING_FORMAT = {
    "record_name" : lambda x : f"{x:<6s}",
    "name"        : lambda x : f" {x:<3s}" if len(x) < 4 else f"{x:4s}",
    "alt"         : lambda x : f"{x:1s}",
    "resn"        : lambda x : f"{x:3s}",
    "chain"       : lambda x : f"{x:1s}",
    "resi"        : lambda x : f"{x:4d}",
    "insertion"   : lambda x : f"{x:1s}",
    "x"           : lambda x : f"{x:8.3f}",
    "y"           : lambda x : f"{x:8.3f}",
    "z"           : lambda x : f"{x:8.3f}",
    "occupancy"   : lambda x : f"{x:6.2f}",
    "b"           : lambda x : f"{x:6.2f}",
    "segi"        : lambda x : f"{x:<3s}",
    "e"           : lambda x : f"{x:2s}",
    "q"           : lambda x : f"{x:2s}"
}