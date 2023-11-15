TEMPERATURE = 300.0          # K
BOLTZMANN   = 1.987204259e-3 # kcal/mol/K (see : 'https://en.wikipedia.org/wiki/Boltzmann_constant')
AVOGADRO    = 6.02214076e23  # 1/mol      (see : 'https://en.wikipedia.org/wiki/Avogadro_constant')

ATOMIC_MASSES = {
    " H" :  1.0079,
    " C" : 12.0107,
    " N" : 14.0067,
    " O" : 15.9994,
    " S" : 32.0650,

    "NA" : 22.9897, # Sodium
    "MG" : 24.3050, # Magnesium
    "CL" : 35.4530, # Chlorine
    " K" : 39.0983, # Potassium
    "AR" : 39.9480, # Argon
    "CA" : 40.0780, # Calcium
    "CU" : 63.5460, # Copper
} # g/mol (see 'https://www.lenntech.com/periodic/mass/atomic-mass.htm')

ATOM_NAME_SELECTION = {
    "ALA" : [" CA "],
    "ARG" : [" N  ", " CA ", " CB ", " CG ", " CD ", " NE ", " CZ ", " NH1"],
    "ASN" : [" N  ", " CA ", " CB ", " CG ", " OD1"],
    "ASP" : [" N  ", " CA ", " CB ", " CG ", " OD1"],
    "CYS" : [" N  ", " CA ", " CB ", " SG "],
    "GLN" : [" N  ", " CA ", " CB ", " CG ", " CD ", " OE1"],
    "GLU" : [" N  ", " CA ", " CB ", " CG ", " CD ", " OE1"],
    "GLY" : [" CA "],
    "HIS" : [" N  ", " CA ", " CB ", " CG ", " ND1"],
    "ILE" : [" N  ", " CA ", " CB ", " CG1", " CD "],
    "LEU" : [" N  ", " CA ", " CB ", " CG ", " CD1"],
    "LYS" : [" N  ", " CA ", " CB ", " CG ", " CD ", " CE ", " NZ "],
    "MET" : [" N  ", " CA ", " CB ", " CG ", " SD ", " CE"],
    "PHE" : [" N  ", " CA ", " CB ", " CG ", " CD1"],
    "PRO" : [" N  ", " CA ", " CB ", " CG ", " CD "],
    "SER" : [" N  ", " CA ", " CB ", " OG "],
    "THR" : [" N  ", " CA ", " CB ", " OG1"],
    "TRP" : [" N  ", " CA ", " CB ", " CG ", " CD1"],
    "TYR" : [" N  ", " CA ", " CB ", " CG ", " CD1"],
    "VAL" : [" N  ", " CA ", " CB ", " CG1"]
} # atoms used for computation of dihedral angles (see 'http://www.mlb.co.jp/linux/science/garlic/doc/commands/dihedrals.html')

from pandas import notnull
PDB_ATOM_FORMAT =  [
    {"id": "record_name", "line": [0, 6]  , "type": str  , "strf": lambda x: "%-6s" % x},
    {"id": "atom_number", "line": [6, 11] , "type": int  , "strf": lambda x: "%+5s" % str(x)},
    {"id": "blank_1"    , "line": [11, 12], "type": str  , "strf": lambda x: "%-1s" % x},
    {"id": "atom_name"  , "line": [12, 16], "type": str  , "strf": lambda x: " %-3s" % x if len(x) < 4 else "%-4s" % x},
    {"id": "alt_loc"    , "line": [16, 17], "type": str  , "strf": lambda x: "%-1s" % x},
    {"id": "resn"       , "line": [17, 20], "type": str  , "strf": lambda x: "%+3s" % x},
    {"id": "blank_2"    , "line": [20, 21], "type": str  , "strf": lambda x: "%-1s" % x},
    {"id": "chain"      , "line": [21, 22], "type": str  , "strf": lambda x: "%-1s" % x},
    {"id": "resi"       , "line": [22, 26], "type": int  , "strf": lambda x: "%+4s" % str(x)},
    {"id": "insertion"  , "line": [26, 27], "type": str  , "strf": lambda x: "%-1s" % x},
    {"id": "blank_3"    , "line": [27, 30], "type": str  , "strf": lambda x: "%-3s" % x},
    {"id": "x"          , "line": [30, 38], "type": float, "strf": lambda x: ("%+8.3f" % x).replace("+", " ")},
    {"id": "y"          , "line": [38, 46], "type": float, "strf": lambda x: ("%+8.3f" % x).replace("+", " ")},
    {"id": "z"          , "line": [46, 54], "type": float, "strf": lambda x: ("%+8.3f" % x).replace("+", " ")},
    {"id": "occupancy"  , "line": [54, 60], "type": float, "strf": lambda x: ("%+6.2f" % x).replace("+", " ")},
    {"id": "b"          , "line": [60, 66], "type": float, "strf": lambda x: ("%+6.2f" % x).replace("+", " ") if len(str(int(x))) < 3 else ("%+6.2f" % x).replace("+", "")},
    {"id": "blank_4"    , "line": [66, 72], "type": str  , "strf": lambda x: "%-7s" % x},
    {"id": "segi"       , "line": [72, 76], "type": str  , "strf": lambda x: "%-3s" % x},
    {"id": "elem"       , "line": [76, 78], "type": str  , "strf": lambda x: "%+2s" % x},
    {"id": "charge"     , "line": [78, 80], "type": float, "strf": lambda x: (("%+2.1f" % x).replace("+", " ") if notnull(x) else "")},
]