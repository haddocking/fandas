import sys

import numpy as np

from fandas.modules.chemical import (
    STANDARD_DATA,
    amino_acids,
    ATOM_LIST,
    c_atm_ind,
    n_atm_ind,
    secondary_structures,
)

import logging

log = logging.getLogger("fandaslog")


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


def replace_bmrb(chem_shifts, bmrb_tables_file, bmrb_columns):
    """Replace the average shifts with the shifts provided in the BMRB table.
    
    Parameters
    ----------
    chem_shifts : numpy.ndarray
        Numpy array containing the chemichal shifts
    bmrb_tables_file : string
        String relates to the PATH where the of the table
    bmrb_columns : list
        List with the indexes (starting-1) of resnum, atom and cshift

    Returns
    -------
    numpy.ndarray
        Numpy array with replaced chemichal shifts

    """
    #  bmrb_columns starts with index 1 so here we need to -1 all positions
    res_index, atom_name_index, cs_index = map(lambda x:x - 1, map(int, bmrb_columns))

    # Read the provided BMRB file
    with open(bmrb_tables_file, "r") as bmrb_file:
        for line in bmrb_file.readlines():

            # Split so that the whitespaces don't matter
            c_shift = line.split()

            # Use the indexes to get:
            #  the provided chemichal shift
            provided_c_shift = float(c_shift[cs_index])
            #  the resnum
            resnum = c_shift[res_index]
            #  the atom
            atom = c_shift[atom_name_index]

            # Find the position of this atom in the chem_shifts table
            position_to_be_replaced = ATOM_LIST.index(atom)

            # Show in the log
            log.info(f"+ {resnum}\t{atom}\t{chem_shifts[res_index][position_to_be_replaced]}\t>> {provided_c_shift}")

            # Replace the appropriate position in the chem_shifts with what we got from the BMRB table
            chem_shifts[res_index][position_to_be_replaced] = provided_c_shift

    return chem_shifts


def assign_chemical_shifts(sequence, sec_struc):
    # reference = np.fromfile("%s/standard.dat" % w_dir, dtype=float, sep=" ").reshape(
    #     ((4, 20, 59))
    # )
    reference = np.fromfile(STANDARD_DATA, dtype=float, sep=" ").reshape(((4, 20, 59)))
    shifts = np.zeros((len(sequence), 59))
    for i, residue in enumerate(sequence):
        for j, amino_acid in enumerate(amino_acids):
            if residue == amino_acid:
                for k, secondary_structure in enumerate(secondary_structures):
                    if sec_struc[i] == secondary_structure:
                        shifts[i] = reference[k, j]
    return shifts


def check_user_input(user_input, input_type, error_type):
    matches = 0
    for component in user_input:
        for reference_component in input_type:
            if component == reference_component:
                matches += 1
    if matches != len(user_input):
        log.error("%s Error: Please check!" % error_type)
        sys.exit()
