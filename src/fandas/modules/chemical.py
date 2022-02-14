from pathlib import Path
import logging

log = logging.getLogger("fandaslog")

STANDARD_DATA = Path(Path(__file__).parents[1], "data/standard.csv")

amino_acids = [
    "A",
    "R",
    "D",
    "N",
    "C",
    "E",
    "Q",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
]

secondary_structures = [
    "a",
    "b",
    "c",
    "n",
]  # a=alpha helix, b=beta sheet, c=random coil, n=bmrb collected average

ATOM_LIST = [
    "H",
    "HA",
    "HA2",
    "HA3",
    "HB",
    "HB2",
    "HB3",
    "HG1",
    "HG2",
    "HG21",
    "HG22",
    "HG3",
    "HD",
    "HD1",
    "HD2",
    "HD21",
    "HD22",
    "HD3",
    "HE",
    "HE1",
    "HE12",
    "HE13",
    "HE2",
    "HE3",
    "HH",
    "HH11",
    "HH12",
    "HH2",
    "HH21",
    "HH22",
    "HZ",
    "HZ2",
    "HZ3",
    "N",
    "NG1",
    "NG2",
    "ND",
    "ND1",
    "ND2",
    "NH1",
    "NH2",
    "NZ",
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

atom_type = [
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "N",
    "N",
    "N",
    "N",
    "N",
    "N",
    "N",
    "N",
    "N",
    "C",
    "C",
    "C",
    "C",
    "C",
    "C",
    "C",
    "C",
    "C",
    "C",
    "C",
    "C",
    "C",
    "C",
    "C",
    "C",
    "C",
]

atom_positions = [
    "N",
    "A",
    "A2",
    "A3",
    "B",
    "B2",
    "B3",
    "G1",
    "G2",
    "G2",
    "G2",
    "G3",
    "D",
    "D1",
    "D2",
    "D2",
    "D2",
    "D3",
    "E",
    "E1",
    "E1",
    "E1",
    "E2",
    "E3",
    "H",
    "H1",
    "H1",
    "H2",
    "H2",
    "H2",
    "Z",
    "Z2",
    "Z3",
    "N",
    "G1",
    "G2",
    "D",
    "D1",
    "D2",
    "H1",
    "H2",
    "Z",
    "C",
    "A",
    "B",
    "G",
    "G1",
    "G2",
    "D",
    "D1",
    "D2",
    "D3",
    "E",
    "E1",
    "E2",
    "H2",
    "Z",
    "Z2",
    "Z3",
]

simple_atom_positions = [
    "N",
    "A",
    "A",
    "A",
    "B",
    "B",
    "B",
    "G",
    "G",
    "G",
    "G",
    "G",
    "D",
    "D",
    "D",
    "D",
    "D",
    "D",
    "E",
    "E",
    "E",
    "E",
    "E",
    "E",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "Z",
    "Z",
    "Z",
    "N",
    "G",
    "G",
    "D",
    "D",
    "D",
    "H",
    "H",
    "Z",
    "C",
    "A",
    "B",
    "G",
    "G",
    "G",
    "D",
    "D",
    "D",
    "D",
    "E",
    "E",
    "E",
    "H",
    "Z",
    "Z",
    "Z",
]

h_atm_ind = [
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
    17,
    18,
    19,
    20,
    21,
    22,
    23,
    24,
    25,
    26,
    27,
    28,
    29,
    30,
    31,
    32,
]

n_atm_ind = [33, 34, 35, 36, 37, 38, 39, 40, 41]

c_atm_ind = [42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58]

list_2d = [
    "NH",
    "HN",
    "CH",
    "HC",
    "HH",
    "DQSQ",
    "CC_SPINDIFF_INTRA",
    "CC_SPINDIFF_INTER",
    "NCA",
    "NCO",
    "NCACX",
    "NCACX_INTER",
    "NCOCX",
    "NCOCA_CB",
    "CANH",
    "CONH",
    "CACONH",
    "COCANH",
    "NCAH",
]

list_3d = [
    "NCACX",
    "NCACX_INTER",
    "NCOCX",
    "NCOCA_CB",
    "SQSQSQ_INTER",
    "DQSQSQ_INTRA",
    "DQSQSQ_INTER",
    "CANH",
    "CONH",
    "CACONH",
    "COCANH",
    "NCAH",
]

list_2dd = ["CC_SPINDIFF", "HH", "CHHC", "NHHC", "CHH", "NHH", "HHC", "NCACX", "NCOCX"]

list_3dd = ["CHH", "NHH", "NCACX", "NCOCX"]
