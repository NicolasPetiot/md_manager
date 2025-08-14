from warnings import filterwarnings
filterwarnings("ignore", category=UserWarning, module="MDAnalysis")

__all__ = [
    "TEMPERATURE",
    "BOLTZMANN",
    "AVOGADRO",
    "ONE_LETTER_CODE",
    "THREE_LETTER_CODE",
    "ATOM_NAME_CHI",
]

TEMPERATURE = 300.0          # K
BOLTZMANN   = 1.987204259e-3 # kcal/mol/K (see : 'https://en.wikipedia.org/wiki/Boltzmann_constant')
AVOGADRO    = 6.02214076e23  # 1/mol      (see : 'https://en.wikipedia.org/wiki/Avogadro_constant')

ONE_LETTER_CODE ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q',
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',
'GLY':'G', 'PRO':'P', 'CYS':'C'}

THREE_LETTER_CODE ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN',
'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR',
'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA',
'G':'GLY', 'P':'PRO', 'C':'CYS'}

ATOM_NAME_CHI = ["N", "CA", "CB", "CG", "SG", "CG1", "OG1", "CD", "SD", "CD1", "OD1", "ND1", "CE", "NE", "OE1", "CZ", "NZ", "NH1"]

ATTRIBUTE_RECORD_EQUIVALENCE = [
    ("record_types", "record_name"),
    ("names", "name"),
    ("altLocs", "alt"),
    ("resnames", "resn"),
    ("chainIDs", "chain"),
    ("resids", "resi"),
    ("icodes", "icode"),
    ("occupancies", "occupancy"),
    ("tempfactors", "b"),
    ("segindices", "segi"),
    ("elements", "e"),
    ("charges", "q"),
    ("masses", "m"),
    ("types", "type"),
    ("ids", "atom_id")
]