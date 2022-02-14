import sys
import copy
import numpy as np
import pandas as pd
from fandas.modules.chemical import (
    STANDARD_DATA,
    ATOM_LIST,
    c_atm_ind,
    n_atm_ind,
)

import logging

log = logging.getLogger("fandaslog")

SS_REF = {"n": "AVG", "a": "ALPHA", "b": "BETA", "c": "COIL"}

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


def write_2d(shift_a, shift_b, extension, sl=True):
    peak = []
    for i in range(len(shift_a)):
        column_2 = round(shift_a[i][2], 4)
        column_3 = round(shift_b[i][2], 4)
        column_4 = "%s%s%s-%s%s%s" % (
            shift_a[i][0],
            shift_a[i][1],
            shift_a[i][3],
            shift_b[i][0],
            shift_b[i][1],
            shift_b[i][3],
        )
        if sl is True:
            column_1 = column_4
        else:
            column_1 = "?-?"
        peak.append([column_1, column_2, column_3, column_4])
    log.info(f"Writing {extension}_exp.txt")
    with open(f"{extension}_exp.txt", "w") as output:
        for pk in peak:
            output.write("%s\t%s\t%s\t%s\n" % (pk[0], pk[1], pk[2], pk[3]))


def write_2dd(shift_a, shift_b, extension, distance_label, sl=True):
    peak = []
    for i in range(len(shift_a)):
        column_2 = round(shift_a[i][2], 4)
        column_3 = round(shift_b[i][2], 4)
        column_4 = "%s%s%s-%s%s%s_%s" % (
            shift_a[i][0],
            shift_a[i][1],
            shift_a[i][3],
            shift_b[i][0],
            shift_b[i][1],
            shift_b[i][3],
            distance_label[i],
        )
        if sl is True:
            column_1 = column_4
        else:
            column_1 = "?-?"
        peak.append([column_1, column_2, column_3, column_4])
    log.info(f"Writing {extension}_exp.txt")
    with open(f"{extension}_exp.txt", "w") as output:
        for pk in peak:
            output.write("%s\t%s\t%s\t%s\n" % (pk[0], pk[1], pk[2], pk[3]))


def write_3dd(shift_a, shift_b, shift_c, extension, distance_label, sl=True):
    peak = []
    for i in range(len(shift_a)):
        column_2 = round(shift_a[i][2], 4)
        column_3 = round(shift_b[i][2], 4)
        column_4 = round(shift_c[i][2], 4)
        column_5 = "%s%s%s-%s%s%s-%s%s%s_%s" % (
            shift_a[i][0],
            shift_a[i][1],
            shift_a[i][3],
            shift_b[i][0],
            shift_b[i][1],
            shift_b[i][3],
            shift_c[i][0],
            shift_c[i][1],
            shift_c[i][3],
            distance_label[i],
        )
        if sl is True:
            column_1 = column_5
        else:
            column_1 = "?-?-?"
        peak.append([column_1, column_2, column_3, column_4, column_5])
    log.info(f"Writing {extension}_exp.txt")
    with open(f"{extension}_exp.txt", "w") as output:
        for pk in peak:
            output.write("%s\t%s\t%s\t%s\t%s\n" % (pk[0], pk[1], pk[2], pk[3], pk[4]))


def write_3d(shift_a, shift_b, shift_c, extension, sl=True):
    peak = []
    for i in range(len(shift_a)):
        column_2 = round(shift_a[i][2], 4)
        column_3 = round(shift_b[i][2], 4)
        column_4 = round(shift_c[i][2], 4)
        column_5 = "%s%s%s-%s%s%s-%s%s%s" % (
            shift_a[i][0],
            shift_a[i][1],
            shift_a[i][3],
            shift_b[i][0],
            shift_b[i][1],
            shift_b[i][3],
            shift_c[i][0],
            shift_c[i][1],
            shift_c[i][3],
        )
        if sl is True:
            column_1 = column_5
        else:
            column_1 = "?-?-?"
        peak.append([column_1, column_2, column_3, column_4, column_5])
    log.info(f"Writing {extension}_exp.txt")
    with open(f"{extension}_exp.txt", "w") as output:
        for pk in peak:
            output.write("%s\t%s\t%s\t%s\t%s\n" % (pk[0], pk[1], pk[2], pk[3], pk[4]))


def fractional_deuteration(sequence, chem_shifts):
    for i, residue in enumerate(sequence):
        if residue == "A":
            rem_atm = ["HA"]
        elif residue == "R":
            rem_atm = ["HA", "HB"]
        elif residue == "D":
            rem_atm = ["HA"]
        elif residue == "N":
            rem_atm = ["HA"]
        elif residue == "C":
            rem_atm = ["HA"]
        elif residue == "E":
            rem_atm = ["HA", "HB"]
        elif residue == "Q":
            rem_atm = ["HA", "HB"]
        elif residue == "G":
            rem_atm = ["HA"]
        elif residue == "H":
            rem_atm = ["HA"]
        elif residue == "I":
            rem_atm = ["HA", "HB", "HG"]
        elif residue == "K":
            rem_atm = ["HA"]
        elif residue == "M":
            rem_atm = ["HA"]
        elif residue == "P":
            rem_atm = ["HA", "HB"]
        elif residue == "L":
            rem_atm = ["HA", "HB", "HG"]
        elif residue == "F":
            rem_atm = ["HA"]
        elif residue == "S":
            rem_atm = ["HA"]
        elif residue == "T":
            rem_atm = ["HA"]
        elif residue == "W":
            rem_atm = ["HA"]
        elif residue == "Y":
            rem_atm = ["HA"]
        elif residue == "V":
            rem_atm = ["HA", "HB"]
        for j, atom in enumerate(ATOM_LIST):
            for element in rem_atm:
                if atom == element:
                    chem_shifts[i][j] = 0.00
    return chem_shifts


def glycerol_label(sequence, chem_shifts, label_type):
    if label_type == 13:
        for i, residue in enumerate(sequence):
            if residue == "A":
                rem_atm = [
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
                ]
            elif residue == "R":
                rem_atm = [
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
                ]
            elif residue == "D":
                rem_atm = [
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
                ]
            elif residue == "N":
                rem_atm = [
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
                ]
            elif residue == "C":
                rem_atm = [
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
                ]
            elif residue == "E":
                rem_atm = [
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
                ]
            elif residue == "Q":
                rem_atm = [
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
                ]
            elif residue == "G":
                rem_atm = [
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
                ]
            elif residue == "H":
                rem_atm = [
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
                ]
            elif residue == "I":
                rem_atm = ["CB", "CD2", "CE1", "CE2", "CE3", "CZ1", "CZ2", "CZ3", "CH"]
            elif residue == "L":
                rem_atm = [
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
                ]
            elif residue == "K":
                rem_atm = ["CG2", "CD2", "CE2", "CE3", "CZ1", "CZ2", "CZ3", "CH"]
            elif residue == "M":
                rem_atm = ["CG2", "CD1", "CD2", "CE2", "CE3", "CZ1", "CZ2", "CZ3", "CH"]
            elif residue == "F":
                rem_atm = ["CA", "CG1", "CG2", "CE2", "CE3", "CZ2", "CZ3", "CH"]
            elif residue == "P":
                rem_atm = [
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
                ]
            elif residue == "S":
                rem_atm = [
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
                ]
            elif residue == "T":
                rem_atm = [
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
                ]
            elif residue == "W":
                rem_atm = ["CA", "CG2", "CD2", "CE1", "CZ1", "CZ3"]
            elif residue == "Y":
                rem_atm = ["CA", "CG1", "CG2", "CE2", "CE3", "CZ2", "CZ3", "CH"]
            elif residue == "V":
                rem_atm = [
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
                ]
            for j, atom in enumerate(ATOM_LIST):
                for element in rem_atm:
                    if atom == element:
                        chem_shifts[i][j] = 0.00
    elif label_type == 2:
        for i, residue in enumerate(sequence):
            if residue == "A":
                rem_atm = [
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
                ]
            elif residue == "R":
                rem_atm = [
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
                ]
            elif residue == "D":
                rem_atm = [
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
                ]
            elif residue == "N":
                rem_atm = [
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
                ]
            elif residue == "C":
                rem_atm = [
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
                ]
            elif residue == "E":
                rem_atm = [
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
                ]
            elif residue == "Q":
                rem_atm = [
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
                ]
            elif residue == "G":
                rem_atm = [
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
                ]
            elif residue == "H":
                rem_atm = [
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
                ]
            elif residue == "I":
                rem_atm = ["CG2", "CD2", "CE1", "CE2", "CE3", "CZ1", "CZ2", "CZ3", "CH"]
            elif residue == "L":
                rem_atm = [
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
                ]
            elif residue == "K":
                rem_atm = ["CG2", "CD2", "CE2", "CE3", "CZ1", "CZ2", "CZ3", "CH"]
            elif residue == "M":
                rem_atm = [
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
                ]
            elif residue == "F":
                rem_atm = [
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
                ]
            elif residue == "P":
                rem_atm = [
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
                ]
            elif residue == "S":
                rem_atm = [
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
                ]
            elif residue == "T":
                rem_atm = [
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
                ]
            elif residue == "W":
                rem_atm = [
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
                ]
            elif residue == "Y":
                rem_atm = [
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
                ]
            elif residue == "V":
                rem_atm = [
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
                ]
            for j, atom in enumerate(ATOM_LIST):
                for element in rem_atm:
                    if atom == element:
                        chem_shifts[i][j] = 0.00
    else:
        log.error("ERROR: %s is an invalid glycerol labelling scheme" % label_type)
        sys.exit()
    return chem_shifts


def rev_label(sequence, chem_shifts, order, amino_acid):
    if order == "rev_cn":
        for i, chem_shift in enumerate(chem_shifts):
            if sequence[i] == amino_acid:
                for n_ind in n_atm_ind:
                    chem_shift[n_ind] = 0.00
        for i, chem_shift in enumerate(chem_shifts):
            if sequence[i] == amino_acid:
                for c_ind in c_atm_ind:
                    chem_shift[c_ind] = 0.00
    if order == "rev_n":
        for i, chem_shift in enumerate(chem_shifts):
            if sequence[i] == amino_acid:
                for n_ind in n_atm_ind:
                    chem_shift[n_ind] = 0.00
    if order == "rev_c":
        for i, chem_shift in enumerate(chem_shifts):
            if sequence[i] == amino_acid:
                for c_ind in c_atm_ind:
                    chem_shift[c_ind] = 0.00
    return chem_shifts


def replace_bmrb(chem_shifts, bmrb_tables_file, bmrb_columns, bt_seq_start):
    """Replace the average shifts with the shifts provided in the BMRB table.

    Parameters
    ----------
    chem_shifts : dict
       Dictionary containing chemichal shifts: {1: {'ATOM': value,}, ...}
    bmrb_tables_file : string
        String related to the PATH where the table file is
    bmrb_columns : list
        List with the indexes (starting-1) of resnum, atom and cshift
    bt_seq_start : int
        Offset to match the BMRB table with sequence

    Returns
    -------
    dict
        Dictionary containing the modficied chemichal shifts: {1: {'ATOM': value,}, ...}

    """
    new_chem_shifts = copy.deepcopy(chem_shifts)
    #  bmrb_columns starts with index 1 so here we need to -1 all positions
    res_index, atom_name_index, cs_index = map(lambda x: x - 1, map(int, bmrb_columns))
    log.info(f"Offset={bt_seq_start}")
    log.info("BMRB_num\tid\tatom\told\t\tnew")

    # Read the provided BMRB file
    with open(bmrb_tables_file, "r") as bmrb_file:
        for line in bmrb_file.readlines():

            # Split so that the whitespaces don't matter
            c_shift = line.split()

            # Use the indexes to get:
            #  the provided chemichal shift
            provided_c_shift = float(c_shift[cs_index])
            #  which position in the chem_shifts, take into account the sequence
            table_index = int(c_shift[res_index]) - bt_seq_start + 1
            # resnum = int(c_shift[res_index]) - 1
            #  the atom
            atom = c_shift[atom_name_index]

            ident = chem_shifts[table_index]["id"]
            # provided_ident = f"{c_shift[res_index]}"
            try:
                default_value = chem_shifts[table_index][atom]
                log.info(
                    f"# ({c_shift[res_index]})\t{ident}\t{atom}\t{default_value}\t>>"
                    f" {provided_c_shift}"
                )
                new_chem_shifts[table_index][atom] = provided_c_shift
            except KeyError:
                log.warning(
                    f"# ({c_shift[res_index]})\t{ident}\t{atom} not found in "
                    "standard table, skipping..."
                )

    return new_chem_shifts


def assign_chemical_shifts(sequence, secondary_structure, data=STANDARD_DATA):
    """Load the standard chemichal shift data and filter by sequence and ss.

    Parameters
    ----------
    sequence : str
        String containing the sequence in one-letter format
    secondary_structure : list
        List containing the secondary structure code
    data : pathlib.Path
        Location of the standard data

    Returns
    -------
    dict
       Dictionary containing chemichal shifts: {1: {'ATOM': value,}, ...}

    """
    log.info(f"Loading standard chemichal shifts from {data}")
    df = pd.read_csv(data)
    shifts = {}
    for i, (res, ss) in enumerate(zip(sequence, secondary_structure), start=1):
        sub_values = df[
            (df["RESNAME"] == AA_REF[res]) & (df["SECONDARY_STRUCTURE"] == SS_REF[ss])
        ]
        shifts[i] = dict(
            (e, v) for e, v in zip(df.columns[2:], sub_values.iloc[0].values[2:])
        )
        shifts[i]["id"] = f"{i}.{res}.{ss}"

    return shifts


def dict2array(chemichal_shifts):
    """Convert the chemichal shift dictionary to numpy.array."""
    shifts = []
    for i in chemichal_shifts:
        shifts.append(
            np.array([chemichal_shifts[i][e] for e in chemichal_shifts[i] if e != "id"])
        )
    return np.array(shifts)


def check_user_input(user_input, input_type, error_type):
    matches = 0
    for component in user_input:
        for reference_component in input_type:
            if component == reference_component:
                matches += 1
    if matches != len(user_input):
        log.error("%s Error: Please check!" % error_type)
        sys.exit()
