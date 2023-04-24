import os
from pathlib import Path

import toml

STANDARD_DATA = Path(Path(__file__).parents[1], "data/standard.csv")

ATOM_LIST = open(STANDARD_DATA).readlines()[0].split(os.linesep)[0].split(",")[2:]

ATOM_REF = {
    "CO": "C",
    "NH": "H",
    "-CA,CX": [atom for atom in ATOM_LIST if "C" in atom and "CA" not in atom],
    "CA,CB": [atom for atom in ATOM_LIST if atom in ["CA", "CB"]],
    "CX": [atom for atom in ATOM_LIST if "C" in atom],
    "HX": [atom for atom in ATOM_LIST if "H" in atom],
}

QUANTUM_RELATIONSHIP = [
    ("C+CA", "C"),
    ("CA+C", "CA"),
    ("CA+CB", "CA"),
    ("CB+CA", "CB"),
    ("CB+CG", "CB"),
    # FIXME: should be CG*, CE*, CD*, etc.
    ("CG+CB", "CG"),
    ("CG+CD1", "CG"),
    ("CG+CD2", "CG"),
    ("CD1+CG", "CD1"),
    ("CD1+CE1", "CD1"),
    ("CD1+CE2", "CD1"),
    ("CD2+CG", "CD2"),
    ("CD2+CE1", "CD2"),
    ("CD2+CE2", "CD2"),
    ("CE1+CD1", "CE1"),
    ("CE1+CD2", "CE1"),
    ("CE1+CZ", "CE1"),
    ("CE2+CD1", "CE2"),
    ("CE2+CD2", "CE2"),
    ("CE2+CZ", "CE2"),
    ("CZ+CE1", "CZ"),
    ("CZ+CE2", "CZ"),
]

EXPERIMENT_CATALOG = toml.load(Path(Path(__file__).parents[1], "data/catalog.toml"))

# This tells the position of each atom in the chemical shift:
#  ex:
# ATOM_INDEXES['HB2']
#  > 5
ATOM_INDEX_DIC = dict((atom, i) for i, atom in enumerate(ATOM_LIST))


N_ATOMS = ["N", "NG1", "NG2", "ND", "ND1", "ND2", "NH1", "NH2", "NZ"]


C_ATOMS = [
    "C",
    "CA",
    "CB",
    "CG",
    "CG1",
    "CG2",
    "CD",
    "CD1",
    "CD2",
    "CD3",
    "CE",
    "CE1",
    "CE2",
    "CH2",
    "CZ",
    "CZ2",
    "CZ3",
]


# Atoms to be removed from each aminoacid based on 1,3,13C,glycerol labeling
#  all atoms present on standard.csv are represented here, even if non existent
#  such as CG1 for ALA
GL_13 = {
    "A": [
        "CA",
        "CG1",
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "R": [
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "D": [
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "N": [
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "C": [
        "CA",
        "CG1",
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "E": [
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "Q": [
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "G": [
        "CA",
        "CB",
        "CG1",
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "H": [
        "CA",
        "CG2",
        "CD1",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "I": ["CB", "CD2", "CE1", "CE2", "CE3", "CZ1", "CZ2", "CZ3", "CH"],
    "L": [
        "C",
        "CB",
        "CG1",
        "CG2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "K": ["CG2", "CD2", "CE2", "CE3", "CZ1", "CZ2", "CZ3", "CH"],
    "M": ["CG2", "CD1", "CD2", "CE2", "CE3", "CZ1", "CZ2", "CZ3", "CH"],
    "F": ["CA", "CG1", "CG2", "CE2", "CE3", "CZ2", "CZ3", "CH"],
    "P": [
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "S": [
        "CA",
        "CG1",
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "T": [
        "CG1",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "W": ["CA", "CG2", "CD2", "CE1", "CZ1", "CZ3"],
    "Y": ["CA", "CG1", "CG2", "CE2", "CE3", "CZ2", "CZ3", "CH"],
    "V": [
        "CA",
        "CB",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
}

AA_REF = {
    "A": "ALA",
    "R": "ARG",
    "D": "ASP",
    "N": "ASN",
    "C": "CYS",
    "E": "GLU",
    "Q": "GLN",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}

SS_REF = {"n": "COIL", "a": "ALPHA", "b": "BETA", "c": "COIL"}

GL_2 = {
    "A": [
        "C",
        "CB",
        "CG1",
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "R": [
        "CG1",
        "CG2",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "D": [
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "N": [
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "C": [
        "C",
        "CB",
        "CG1",
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "E": [
        "CG1",
        "CG2",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "Q": [
        "CG1",
        "CG2",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "G": [
        "C",
        "CB",
        "CG1",
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "H": [
        "C",
        "CB",
        "CG2",
        "CD1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "I": ["CG2", "CD2", "CE1", "CE2", "CE3", "CZ1", "CZ2", "CZ3", "CH"],
    "L": [
        "CA",
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "K": ["CG2", "CD2", "CE2", "CE3", "CZ1", "CZ2", "CZ3", "CH"],
    "M": [
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "F": [
        "C",
        "CB",
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "P": [
        "CG1",
        "CG2",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "S": [
        "C",
        "CB",
        "CG1",
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "T": [
        "CG1",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "W": [
        "C",
        "CB",
        "CG2",
        "CD1",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CH",
    ],
    "Y": [
        "C",
        "CB",
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
    "V": [
        "C",
        "CG1",
        "CG2",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CE3",
        "CZ1",
        "CZ2",
        "CZ3",
        "CH",
    ],
}

DEUTERATION = {
    "A": ["HA"],
    "R": ["HA", "HB"],
    "D": ["HA"],
    "N": ["HA"],
    "C": ["HA"],
    "E": ["HA", "HB"],
    "Q": ["HA", "HB"],
    "G": ["HA"],
    "H": ["HA"],
    "I": ["HA", "HB", "HG"],
    "K": ["HA"],
    "M": ["HA"],
    "P": ["HA", "HB"],
    "L": ["HA", "HB", "HG"],
    "F": ["HA"],
    "S": ["HA"],
    "T": ["HA"],
    "W": ["HA"],
    "Y": ["HA"],
    "V": ["HA", "HB"],
}


SPECIAL_CASES = [
    "C-C DQ-SQ Correlation",
    "SQSQSQ (residues i, i+1 & i-1)",
    "DQSQSQ intra residue",
    "DQSQSQ (residues i, i+1 & i-1)",
    "C-H",
]
