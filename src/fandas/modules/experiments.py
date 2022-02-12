import logging

from fandas.modules.chemical import (
    atom_positions,
    atom_type,
    ATOM_LIST,
    c_atm_ind,
    h_atm_ind,
    n_atm_ind,
    simple_atom_positions,
)
from fandas.modules.utils import write_2d, write_2dd, write_3d, write_3dd

log = logging.getLogger("fandaslog")


def hh(sequence, chem_shifts):
    # this function produces peak list for intra-residue HH experiment
    shift_1 = []
    shift_2 = []
    for i, residue in enumerate(sequence):  # loop through the sequence
        for j, atom in enumerate(atom_type):  # loop through all atomtypes
            if atom == "H":  # check if it is a carbon atom
                for h_ind in h_atm_ind:  # loop through all the carbon indices
                    if (
                        h_ind != j
                    ):  # ensure that the two indices "j" and "h_ind" are not the same
                        if (
                            chem_shifts[i][h_ind] != 0 and chem_shifts[i][j] != 0
                        ):  # check if the chemical shifts of the two are not zero
                            shift_1.append(
                                [residue, i + 1, chem_shifts[i][j], ATOM_LIST[j]]
                            )
                            shift_2.append(
                                [residue, i + 1, chem_shifts[i][h_ind], ATOM_LIST[h_ind]]
                            )
    write_2d(shift_1, shift_2, "hh")


def hh_dist(sequence, chem_shifts, dlist, distances):
    shift_1 = []
    shift_2 = []
    distance_label = []
    for dist_line in dlist:
        resi_1 = int(dist_line[0]) - 1
        atm_1 = dist_line[1]
        resi_2 = int(dist_line[2]) - 1
        atm_2 = dist_line[3]
        pos_1 = ATOM_LIST.index(atm_1)
        pos_2 = ATOM_LIST.index(atm_2)
        if (
            (atom_type[pos_1] == "H")
            and (atom_type[pos_2] == "H")
            and (chem_shifts[resi_1][pos_1] != 0)
            and (chem_shifts[resi_2][pos_2] != 0)
        ):
            shift_1.append(
                [
                    sequence[resi_1],
                    (resi_1) + 1,
                    chem_shifts[resi_1][pos_1],
                    ATOM_LIST[pos_1],
                ]
            )
            shift_2.append(
                [
                    sequence[resi_2],
                    (resi_2) + 1,
                    chem_shifts[resi_2][pos_2],
                    ATOM_LIST[pos_2],
                ]
            )
            shift_1.append(
                [
                    sequence[resi_2],
                    (resi_2) + 1,
                    chem_shifts[resi_2][pos_2],
                    ATOM_LIST[pos_2],
                ]
            )
            shift_2.append(
                [
                    sequence[resi_1],
                    (resi_1) + 1,
                    chem_shifts[resi_1][pos_1],
                    ATOM_LIST[pos_1],
                ]
            )
            distance_label.append(distances[dlist.index(dist_line)])
            distance_label.append(distances[dlist.index(dist_line)])
    write_2dd(shift_1, shift_2, "hh_dist", distance_label)


def chh_dist_3d(sequence, chem_shifts, dlist, distances):
    shift_1 = []
    shift_2 = []
    shift_3 = []
    # shift_4 = []
    distance_label = []
    for dist_line in dlist:
        resi_1 = int(dist_line[0]) - 1
        atm_1 = dist_line[1]
        resi_2 = int(dist_line[2]) - 1
        atm_2 = dist_line[3]
        pos_1 = ATOM_LIST.index(atm_1)
        pos_2 = ATOM_LIST.index(atm_2)
        if (
            (atom_type[pos_1] == "H")
            and (atom_type[pos_2] == "H")
            and (chem_shifts[resi_1][pos_1] != 0)
            and (chem_shifts[resi_2][pos_2] != 0)
        ):
            for c_ind_1 in c_atm_ind:
                if (
                    simple_atom_positions[c_ind_1] == simple_atom_positions[pos_1]
                ) and (chem_shifts[resi_1][c_ind_1] != 0):
                    for c_ind_2 in c_atm_ind:
                        if (
                            simple_atom_positions[c_ind_2]
                            == simple_atom_positions[pos_2]
                        ) and (chem_shifts[resi_1][c_ind_2] != 0):
                            shift_1.append(
                                [
                                    sequence[resi_1],
                                    (resi_1) + 1,
                                    chem_shifts[resi_1][c_ind_1],
                                    ATOM_LIST[c_ind_1],
                                ]
                            )
                            shift_2.append(
                                [
                                    sequence[resi_1],
                                    (resi_1) + 1,
                                    chem_shifts[resi_1][pos_1],
                                    ATOM_LIST[pos_1],
                                ]
                            )
                            shift_3.append(
                                [
                                    sequence[resi_2],
                                    (resi_2) + 1,
                                    chem_shifts[resi_2][pos_2],
                                    ATOM_LIST[pos_2],
                                ]
                            )
                            shift_1.append(
                                [
                                    sequence[resi_2],
                                    (resi_2) + 1,
                                    chem_shifts[resi_2][c_ind_2],
                                    ATOM_LIST[c_ind_2],
                                ]
                            )
                            shift_2.append(
                                [
                                    sequence[resi_2],
                                    (resi_2) + 1,
                                    chem_shifts[resi_2][pos_2],
                                    ATOM_LIST[pos_2],
                                ]
                            )
                            shift_3.append(
                                [
                                    sequence[resi_1],
                                    (resi_1) + 1,
                                    chem_shifts[resi_1][pos_1],
                                    ATOM_LIST[pos_1],
                                ]
                            )
                            distance_label.append(distances[dlist.index(dist_line)])
                            distance_label.append(distances[dlist.index(dist_line)])
                            distance_label.append(distances[dlist.index(dist_line)])
    write_3dd(shift_1, shift_2, shift_3, "chh_dist_3d", distance_label)


def nhh_dist_3d(sequence, chem_shifts, dlist, distances):
    shift_1 = []
    shift_2 = []
    shift_3 = []
    # shift_4 = []
    distance_label = []
    for dist_line in dlist:
        resi_1 = int(dist_line[0]) - 1
        atm_1 = dist_line[1]
        resi_2 = int(dist_line[2]) - 1
        atm_2 = dist_line[3]
        pos_1 = ATOM_LIST.index(atm_1)
        pos_2 = ATOM_LIST.index(atm_2)
        if (
            (atom_type[pos_1] == "H")
            and (atom_type[pos_2] == "H")
            and (chem_shifts[resi_1][pos_1] != 0)
            and (chem_shifts[resi_2][pos_2] != 0)
        ):
            for n_ind_1 in n_atm_ind:
                if (
                    simple_atom_positions[n_ind_1] == simple_atom_positions[pos_1]
                ) and (chem_shifts[resi_1][n_ind_1] != 0):
                    for n_ind_2 in n_atm_ind:
                        if (
                            simple_atom_positions[n_ind_2]
                            == simple_atom_positions[pos_2]
                        ) and (chem_shifts[resi_1][n_ind_2] != 0):
                            shift_1.append(
                                [
                                    sequence[resi_1],
                                    (resi_1) + 1,
                                    chem_shifts[resi_1][n_ind_1],
                                    ATOM_LIST[n_ind_1],
                                ]
                            )
                            shift_2.append(
                                [
                                    sequence[resi_1],
                                    (resi_1) + 1,
                                    chem_shifts[resi_1][pos_1],
                                    ATOM_LIST[pos_1],
                                ]
                            )
                            shift_3.append(
                                [
                                    sequence[resi_2],
                                    (resi_2) + 1,
                                    chem_shifts[resi_2][pos_2],
                                    ATOM_LIST[pos_2],
                                ]
                            )
                            shift_1.append(
                                [
                                    sequence[resi_2],
                                    (resi_2) + 1,
                                    chem_shifts[resi_2][n_ind_2],
                                    ATOM_LIST[n_ind_2],
                                ]
                            )
                            shift_2.append(
                                [
                                    sequence[resi_2],
                                    (resi_2) + 1,
                                    chem_shifts[resi_2][pos_2],
                                    ATOM_LIST[pos_2],
                                ]
                            )
                            shift_3.append(
                                [
                                    sequence[resi_1],
                                    (resi_1) + 1,
                                    chem_shifts[resi_1][pos_1],
                                    ATOM_LIST[pos_1],
                                ]
                            )
                            distance_label.append(distances[dlist.index(dist_line)])
                            distance_label.append(distances[dlist.index(dist_line)])
                            distance_label.append(distances[dlist.index(dist_line)])
    write_3dd(shift_1, shift_2, shift_3, "nhh_dist_3d", distance_label)


def nhh_dist(sequence, chem_shifts, dlist, distances):
    shift_1 = []
    shift_2 = []
    distance_label = []
    for dist_line in dlist:
        resi_1 = int(dist_line[0]) - 1
        atm_1 = dist_line[1]
        resi_2 = int(dist_line[2]) - 1
        atm_2 = dist_line[3]
        pos_1 = ATOM_LIST.index(atm_1)
        pos_2 = ATOM_LIST.index(atm_2)
        if (
            (atom_type[pos_1] == "H")
            and (atom_type[pos_2] == "H")
            and (chem_shifts[resi_1][pos_1] != 0)
            and (chem_shifts[resi_2][pos_2] != 0)
        ):
            for n_ind_1 in n_atm_ind:
                if (
                    simple_atom_positions[n_ind_1] == simple_atom_positions[pos_1]
                ) and (chem_shifts[resi_1][n_ind_1] != 0):
                    for h_ind_2 in h_atm_ind:
                        if (
                            simple_atom_positions[h_ind_2]
                            == simple_atom_positions[pos_2]
                        ) and (chem_shifts[resi_2][h_ind_2] != 0):
                            shift_1.append(
                                [
                                    sequence[resi_1],
                                    (resi_1) + 1,
                                    chem_shifts[resi_1][n_ind_1],
                                    ATOM_LIST[n_ind_1],
                                ]
                            )
                            shift_2.append(
                                [
                                    sequence[resi_2],
                                    (resi_2) + 1,
                                    chem_shifts[resi_2][h_ind_2],
                                    ATOM_LIST[h_ind_2],
                                ]
                            )
                            shift_1.append(
                                [
                                    sequence[resi_2],
                                    (resi_2) + 1,
                                    chem_shifts[resi_2][h_ind_2],
                                    ATOM_LIST[h_ind_2],
                                ]
                            )
                            shift_2.append(
                                [
                                    sequence[resi_1],
                                    (resi_1) + 1,
                                    chem_shifts[resi_1][n_ind_1],
                                    ATOM_LIST[n_ind_1],
                                ]
                            )
                            distance_label.append(distances[dlist.index(dist_line)])
                            distance_label.append(distances[dlist.index(dist_line)])
    write_2dd(shift_1, shift_2, "nhh_dist_2d", distance_label)


def nhhc_dist(sequence, chem_shifts, dlist, distances):
    shift_1 = []
    shift_2 = []
    distance_label = []
    for dist_line in dlist:
        resi_1 = int(dist_line[0]) - 1
        atm_1 = dist_line[1]
        resi_2 = int(dist_line[2]) - 1
        atm_2 = dist_line[3]
        pos_1 = ATOM_LIST.index(atm_1)
        pos_2 = ATOM_LIST.index(atm_2)
        if (
            (atom_type[pos_1] == "H")
            and (atom_type[pos_2] == "H")
            and (chem_shifts[resi_1][pos_1] != 0)
            and (chem_shifts[resi_2][pos_2] != 0)
        ):
            for n_ind_1 in n_atm_ind:
                if (
                    simple_atom_positions[n_ind_1] == simple_atom_positions[pos_1]
                ) and (chem_shifts[resi_1][n_ind_1] != 0):
                    for c_ind_2 in c_atm_ind:
                        if (
                            simple_atom_positions[c_ind_2]
                            == simple_atom_positions[pos_2]
                        ) and (chem_shifts[resi_2][c_ind_2] != 0):
                            shift_1.append(
                                [
                                    sequence[resi_1],
                                    (resi_1) + 1,
                                    chem_shifts[resi_1][n_ind_1],
                                    ATOM_LIST[n_ind_1],
                                ]
                            )
                            shift_2.append(
                                [
                                    sequence[resi_2],
                                    (resi_2) + 1,
                                    chem_shifts[resi_2][c_ind_2],
                                    ATOM_LIST[c_ind_2],
                                ]
                            )
                            shift_1.append(
                                [
                                    sequence[resi_2],
                                    (resi_2) + 1,
                                    chem_shifts[resi_2][c_ind_2],
                                    ATOM_LIST[c_ind_2],
                                ]
                            )
                            shift_2.append(
                                [
                                    sequence[resi_1],
                                    (resi_1) + 1,
                                    chem_shifts[resi_1][n_ind_1],
                                    ATOM_LIST[n_ind_1],
                                ]
                            )
                            distance_label.append(distances[dlist.index(dist_line)])
                            distance_label.append(distances[dlist.index(dist_line)])
    write_2dd(shift_1, shift_2, "nhhc_dist", distance_label)


def chh_dist(sequence, chem_shifts, dlist, distances):
    shift_1 = []
    shift_2 = []
    distance_label = []
    for dist_line in dlist:
        resi_1 = int(dist_line[0]) - 1
        atm_1 = dist_line[1]
        resi_2 = int(dist_line[2]) - 1
        atm_2 = dist_line[3]
        pos_1 = ATOM_LIST.index(atm_1)
        pos_2 = ATOM_LIST.index(atm_2)
        if (
            (atom_type[pos_1] == "H")
            and (atom_type[pos_2] == "H")
            and (chem_shifts[resi_1][pos_1] != 0)
            and (chem_shifts[resi_2][pos_2] != 0)
        ):
            for c_ind_1 in c_atm_ind:
                if (
                    simple_atom_positions[c_ind_1] == simple_atom_positions[pos_1]
                ) and (chem_shifts[resi_1][c_ind_1] != 0):
                    for h_ind_2 in h_atm_ind:
                        if (
                            simple_atom_positions[h_ind_2]
                            == simple_atom_positions[pos_2]
                        ) and (chem_shifts[resi_2][h_ind_2] != 0):
                            shift_1.append(
                                [
                                    sequence[resi_1],
                                    (resi_1) + 1,
                                    chem_shifts[resi_1][c_ind_1],
                                    ATOM_LIST[c_ind_1],
                                ]
                            )
                            shift_2.append(
                                [
                                    sequence[resi_2],
                                    (resi_2) + 1,
                                    chem_shifts[resi_2][h_ind_2],
                                    ATOM_LIST[h_ind_2],
                                ]
                            )
                            shift_1.append(
                                [
                                    sequence[resi_2],
                                    (resi_2) + 1,
                                    chem_shifts[resi_2][h_ind_2],
                                    ATOM_LIST[h_ind_2],
                                ]
                            )
                            shift_2.append(
                                [
                                    sequence[resi_1],
                                    (resi_1) + 1,
                                    chem_shifts[resi_1][c_ind_1],
                                    ATOM_LIST[c_ind_1],
                                ]
                            )
                            distance_label.append(distances[dlist.index(dist_line)])
                            distance_label.append(distances[dlist.index(dist_line)])
    write_2dd(shift_1, shift_2, "chh_dist_2d", distance_label)


def chhc_dist(sequence, chem_shifts, dlist, distances):
    shift_1 = []
    shift_2 = []
    distance_label = []
    for dist_line in dlist:
        resi_1 = int(dist_line[0]) - 1
        atm_1 = dist_line[1]
        resi_2 = int(dist_line[2]) - 1
        atm_2 = dist_line[3]
        pos_1 = ATOM_LIST.index(atm_1)
        pos_2 = ATOM_LIST.index(atm_2)
        if (
            (atom_type[pos_1] == "H")
            and (atom_type[pos_2] == "H")
            and (chem_shifts[resi_1][pos_1] != 0)
            and (chem_shifts[resi_2][pos_2] != 0)
        ):
            for c_ind_1 in c_atm_ind:
                if (
                    simple_atom_positions[c_ind_1] == simple_atom_positions[pos_1]
                ) and (chem_shifts[resi_1][c_ind_1] != 0):
                    for c_ind_2 in c_atm_ind:
                        if (
                            simple_atom_positions[c_ind_2]
                            == simple_atom_positions[pos_2]
                        ) and (chem_shifts[resi_2][c_ind_2] != 0):
                            shift_1.append(
                                [
                                    sequence[resi_1],
                                    (resi_1) + 1,
                                    chem_shifts[resi_1][c_ind_1],
                                    ATOM_LIST[c_ind_1],
                                ]
                            )
                            shift_2.append(
                                [
                                    sequence[resi_2],
                                    (resi_2) + 1,
                                    chem_shifts[resi_2][c_ind_2],
                                    ATOM_LIST[c_ind_2],
                                ]
                            )
                            shift_1.append(
                                [
                                    sequence[resi_2],
                                    (resi_2) + 1,
                                    chem_shifts[resi_2][c_ind_2],
                                    ATOM_LIST[c_ind_2],
                                ]
                            )
                            shift_2.append(
                                [
                                    sequence[resi_1],
                                    (resi_1) + 1,
                                    chem_shifts[resi_1][c_ind_1],
                                    ATOM_LIST[c_ind_1],
                                ]
                            )
                            distance_label.append(distances[dlist.index(dist_line)])
                            distance_label.append(distances[dlist.index(dist_line)])
    write_2dd(shift_1, shift_2, "chhc_dist", distance_label)


def cc_spindiff_intra(sequence, chem_shifts):
    # this function produces peak list for intra-residue CC spin diffusion experiment
    shift_1 = []
    shift_2 = []
    for i, residue in enumerate(sequence):  # loop through the sequence
        for j, atom in enumerate(atom_type):  # loop through all atomtypes
            if atom == "C":  # check if it is a carbon atom
                for c_ind in c_atm_ind:  # loop through all the carbon indices
                    if (
                        (c_ind != j)
                        and (chem_shifts[i][c_ind] != 0)
                        and (chem_shifts[i][j]) != 0
                    ):
                        shift_1.append([residue, i + 1, chem_shifts[i][j], ATOM_LIST[j]])
                        shift_2.append(
                            [residue, i + 1, chem_shifts[i][c_ind], ATOM_LIST[c_ind]]
                        )
    write_2d(shift_1, shift_2, "cc_spindiff_intra")


def cc_spindiff_inter(sequence, chem_shifts):
    shift_1 = []
    shift_2 = []
    for i, residue in enumerate(sequence):
        for j, atom in enumerate(atom_type):
            if atom == "C":
                for c_ind in c_atm_ind:
                    if (
                        (c_ind != j)
                        and (chem_shifts[i][c_ind] != 0)
                        and (chem_shifts[i][j] != 0)
                    ):
                        shift_1.append([residue, i + 1, chem_shifts[i][j], ATOM_LIST[j]])
                        shift_2.append(
                            [residue, i + 1, chem_shifts[i][c_ind], ATOM_LIST[c_ind]]
                        )
                if i != 0:
                    for c_ind in c_atm_ind:
                        if (chem_shifts[i - 1][c_ind] != 0) and (
                            chem_shifts[i][j] != 0
                        ):
                            shift_1.append(
                                [residue, i + 1, chem_shifts[i][j], ATOM_LIST[j]]
                            )
                            shift_2.append(
                                [
                                    sequence[i - 1],
                                    i,
                                    chem_shifts[i - 1][c_ind],
                                    ATOM_LIST[c_ind],
                                ]
                            )
                if i + 1 != len(sequence):
                    for c_ind in c_atm_ind:
                        if chem_shifts[i + 1][c_ind] != 0 and chem_shifts[i][j] != 0:
                            shift_1.append(
                                [residue, i + 1, chem_shifts[i][j], ATOM_LIST[j]]
                            )
                            shift_2.append(
                                [
                                    sequence[i + 1],
                                    i + 2,
                                    chem_shifts[i + 1][c_ind],
                                    ATOM_LIST[c_ind],
                                ]
                            )
    write_2d(shift_1, shift_2, "cc_spindiff_inter")


def cc_spindiff_dist(sequence, chem_shifts, dlist, distances):
    shift_1 = []
    shift_2 = []
    distance_label = []
    for dist_line in dlist:
        resi_1 = int(dist_line[0]) - 1
        atm_1 = dist_line[1]
        resi_2 = int(dist_line[2]) - 1
        atm_2 = dist_line[3]
        pos_1 = ATOM_LIST.index(atm_1)
        pos_2 = ATOM_LIST.index(atm_2)
        if (
            (atom_type[pos_1] == "C")
            and (atom_type[pos_2] == "C")
            and (chem_shifts[resi_1][pos_1] != 0)
            and (chem_shifts[resi_2][pos_2] != 0)
        ):
            shift_1.append(
                [
                    sequence[resi_1],
                    (resi_1) + 1,
                    chem_shifts[resi_1][pos_1],
                    ATOM_LIST[pos_1],
                ]
            )
            shift_2.append(
                [
                    sequence[resi_2],
                    (resi_2) + 1,
                    chem_shifts[resi_2][pos_2],
                    ATOM_LIST[pos_2],
                ]
            )
            shift_1.append(
                [
                    sequence[resi_2],
                    (resi_2) + 1,
                    chem_shifts[resi_2][pos_2],
                    ATOM_LIST[pos_2],
                ]
            )
            shift_2.append(
                [
                    sequence[resi_1],
                    (resi_1) + 1,
                    chem_shifts[resi_1][pos_1],
                    ATOM_LIST[pos_1],
                ]
            )
            distance_label.append(distances[dlist.index(dist_line)])
            distance_label.append(distances[dlist.index(dist_line)])
    write_2dd(shift_1, shift_2, "cc_spindiff_dist", distance_label)


def nca(sequence, chem_shifts):
    # this function produces peak list for intra-residue NCA experiment
    shift_1 = []
    shift_2 = []
    for i, residue in enumerate(sequence):  # loop through the sequence
        if (chem_shifts[i][ATOM_LIST.index("N")] != 0) and (
            chem_shifts[i][ATOM_LIST.index("CA")] != 0
        ):  # check if the chemical shifts of the two are not zero
            shift_1.append([residue, i + 1, chem_shifts[i][ATOM_LIST.index("N")], "N"])
            shift_2.append([residue, i + 1, chem_shifts[i][ATOM_LIST.index("CA")], "CA"])
    write_2d(shift_1, shift_2, "nca")


def nco(sequence, chem_shifts):
    # this function produces peak list for inter-residue NCO experiment
    shift_1 = []
    shift_2 = []
    for i, residue in enumerate(sequence):  # loop through the sequence
        if (
            (i != 0)
            and (chem_shifts[i][ATOM_LIST.index("N")] != 0)
            and (chem_shifts[i - 1][ATOM_LIST.index("C")] != 0)
        ):  # check if the chemical shifts of the two are not zero
            shift_1.append([residue, i + 1, chem_shifts[i][ATOM_LIST.index("N")], "N"])
            shift_2.append(
                [sequence[i - 1], i, chem_shifts[i - 1][ATOM_LIST.index("C")], "C"]
            )
    write_2d(shift_1, shift_2, "nco")


def ncocx_dist(sequence, chem_shifts, dimensionality, dlist, distances):
    shift_1 = []
    shift_2 = []
    shift_3 = []
    distance_label = []
    for dist_line in dlist:
        resi_1 = int(dist_line[0]) - 1
        atm_1 = dist_line[1]
        resi_2 = int(dist_line[2]) - 1
        atm_2 = dist_line[3]
        pos_1 = ATOM_LIST.index(atm_1)
        pos_2 = ATOM_LIST.index(atm_2)
        if (
            (atm_1 == "C")
            and (atom_type[pos_2] == "C")
            and (chem_shifts[resi_1 + 1][ATOM_LIST.index("N")] != 0)
            and (chem_shifts[resi_1][pos_1] != 0)
            and (chem_shifts[resi_2][pos_2] != 0)
        ):
            shift_1.append(
                [
                    sequence[(resi_1) + 1],
                    (resi_1) + 2,
                    chem_shifts[(resi_1) + 1][ATOM_LIST.index("N")],
                    "N",
                ]
            )
            shift_2.append(
                [sequence[resi_1], resi_1 + 1, chem_shifts[resi_1][pos_1], "CA"]
            )
            shift_3.append(
                [sequence[resi_2], resi_2 + 1, chem_shifts[resi_2][pos_2], ATOM_LIST[pos_2]]
            )
            distance_label.append(distances[dlist.index(dist_line)])
        elif (
            (atm_2 == "C")
            and (atom_type[pos_1] == "C")
            and (chem_shifts[resi_2 + 1][ATOM_LIST.index("N")] != 0)
            and (chem_shifts[resi_2][pos_2] != 0)
            and (chem_shifts[resi_1][pos_1] != 0)
        ):
            shift_1.append(
                [
                    sequence[(resi_2) + 1],
                    (resi_2) + 2,
                    chem_shifts[(resi_2) + 1][ATOM_LIST.index("N")],
                    "N",
                ]
            )
            shift_2.append(
                [sequence[resi_2], resi_2 + 1, chem_shifts[resi_2][pos_2], "CA"]
            )
            shift_3.append(
                [sequence[resi_1], resi_1 + 1, chem_shifts[resi_1][pos_1], ATOM_LIST[pos_1]]
            )
            distance_label.append(distances[dlist.index(dist_line)])
    if dimensionality == 2:
        write_2dd(shift_1, shift_3, "ncocx_dist_2d", distance_label)
    elif dimensionality == 3:
        write_3dd(shift_1, shift_2, shift_3, "ncocx_dist_3d", distance_label)


def ncacx_dist(sequence, chem_shifts, dimensionality, dlist, distances):
    shift_1 = []
    shift_2 = []
    shift_3 = []
    distance_label = []
    for dist_line in dlist:
        resi_1 = int(dist_line[0]) - 1
        atm_1 = dist_line[1]
        resi_2 = int(dist_line[2]) - 1
        atm_2 = dist_line[3]
        pos_1 = ATOM_LIST.index(atm_1)
        pos_2 = ATOM_LIST.index(atm_2)
        if (
            (atm_1 == "CA")
            and (chem_shifts[resi_1][ATOM_LIST.index("N")] != 0)
            and (chem_shifts[resi_1][pos_1] != 0)
            and (chem_shifts[resi_2][pos_2] != 0)
        ):
            shift_1.append(
                [
                    sequence[resi_1],
                    resi_1 + 1,
                    chem_shifts[resi_1][ATOM_LIST.index("N")],
                    "N",
                ]
            )
            shift_2.append(
                [sequence[resi_1], resi_1 + 1, chem_shifts[resi_1][pos_1], "CA"]
            )
            shift_3.append(
                [sequence[resi_2], resi_2 + 1, chem_shifts[resi_2][pos_2], ATOM_LIST[pos_2]]
            )
            distance_label.append(distances[dlist.index(dist_line)])
        elif (
            (atm_2 == "CA")
            and (chem_shifts[resi_2][ATOM_LIST.index("N")] != 0)
            and (chem_shifts[resi_2][pos_2] != 0)
            and (chem_shifts[resi_1][pos_1] != 0)
        ):
            shift_1.append(
                [
                    sequence[resi_2],
                    resi_2 + 1,
                    chem_shifts[resi_2][ATOM_LIST.index("N")],
                    "N",
                ]
            )
            shift_2.append(
                [sequence[resi_2], resi_2 + 1, chem_shifts[resi_2][pos_2], "CA"]
            )
            shift_3.append(
                [sequence[resi_1], resi_1 + 1, chem_shifts[resi_1][pos_1], ATOM_LIST[pos_1]]
            )
            distance_label.append(distances[dlist.index(dist_line)])
    if dimensionality == 2:
        write_2dd(shift_1, shift_3, "ncacx_dist_2d", distance_label)
    elif dimensionality == 3:
        write_3dd(shift_1, shift_2, shift_3, "ncacx_dist_3d", distance_label)


def ncacx_inter(sequence, chem_shifts, dimensionality):
    # this function produces peak list for inter residue NCACX experiment
    shift_1 = []
    shift_2 = []
    shift_3 = []
    for i, residue in enumerate(sequence):  # loop through the sequence
        if (chem_shifts[i][ATOM_LIST.index("N")] != 0) and (
            chem_shifts[i][ATOM_LIST.index("CA")] != 0
        ):  # check if the chemical shifts of the two are not zero
            for c_ind in c_atm_ind:
                if (chem_shifts[i][c_ind] != 0) and (ATOM_LIST[c_ind] != "CA"):
                    shift_1.append(
                        [residue, i + 1, chem_shifts[i][ATOM_LIST.index("N")], "N"]
                    )
                    shift_2.append(
                        [residue, i + 1, chem_shifts[i][ATOM_LIST.index("CA")], "CA"]
                    )
                    shift_3.append(
                        [residue, i + 1, chem_shifts[i][c_ind], ATOM_LIST[c_ind]]
                    )
            if i != 0:
                for c_ind in c_atm_ind:
                    if (chem_shifts[i - 1][c_ind] != 0) and (ATOM_LIST[c_ind] != "CA"):
                        shift_1.append(
                            [residue, i + 1, chem_shifts[i][ATOM_LIST.index("N")], "N"]
                        )
                        shift_2.append(
                            [residue, i + 1, chem_shifts[i][ATOM_LIST.index("CA")], "CA"]
                        )
                        shift_3.append(
                            [
                                sequence[i - 1],
                                i,
                                chem_shifts[i - 1][c_ind],
                                ATOM_LIST[c_ind],
                            ]
                        )
            if i + 1 != len(sequence):
                for c_ind in c_atm_ind:
                    if (chem_shifts[i + 1][c_ind] != 0) and (ATOM_LIST[c_ind] != "CA"):
                        shift_1.append(
                            [residue, i + 1, chem_shifts[i][ATOM_LIST.index("N")], "N"]
                        )
                        shift_2.append(
                            [residue, i + 1, chem_shifts[i][ATOM_LIST.index("CA")], "CA"]
                        )
                        shift_3.append(
                            [
                                sequence[i + 1],
                                i + 2,
                                chem_shifts[i + 1][c_ind],
                                ATOM_LIST[c_ind],
                            ]
                        )
    if dimensionality == 2:
        write_2d(shift_1, shift_3, "ncacx_inter_2d")
    elif dimensionality == 3:
        write_3d(shift_1, shift_2, shift_3, "ncacx_inter_3d")


def ncacx(sequence, chem_shifts, dimensionality):
    # this function produces peak list for inter residue NCACX experiment
    shift_1 = []
    shift_2 = []
    shift_3 = []
    for i, residue in enumerate(sequence):  # loop through the sequence
        if (
            (chem_shifts[i][ATOM_LIST.index("H")] != 0)
            and (chem_shifts[i][ATOM_LIST.index("N")] != 0)
            and (chem_shifts[i][ATOM_LIST.index("CA")] != 0)
        ):  # check if the chemical shifts of the two are not zero
            for c_ind in c_atm_ind:
                if (chem_shifts[i][c_ind] != 0) and (ATOM_LIST[c_ind] != "CA"):
                    shift_1.append(
                        [residue, i + 1, chem_shifts[i][ATOM_LIST.index("N")], "N"]
                    )
                    shift_2.append(
                        [residue, i + 1, chem_shifts[i][ATOM_LIST.index("CA")], "CA"]
                    )
                    shift_3.append(
                        [residue, i + 1, chem_shifts[i][c_ind], ATOM_LIST[c_ind]]
                    )
    if dimensionality == 2:
        write_2d(shift_1, shift_3, "ncacx_2d")
    elif dimensionality == 3:
        write_3d(shift_1, shift_2, shift_3, "ncacx_3d")


def ncocx(sequence, chem_shifts, dimensionality):
    # this function produces peak list for inter residue NCOCX experiment
    shift_1 = []
    shift_2 = []
    shift_3 = []
    for i, residue in enumerate(sequence):  # loop through the sequence
        if (
            (i != 0)
            and (chem_shifts[i][ATOM_LIST.index("N")] != 0)
            and (chem_shifts[i - 1][ATOM_LIST.index("C")] != 0)
        ):  # check if the chemical shifts of the two are not zero
            for c_ind in c_atm_ind:
                if (chem_shifts[i - 1][c_ind] != 0) and (ATOM_LIST[c_ind] != "C"):
                    shift_1.append(
                        [residue, i + 1, chem_shifts[i][ATOM_LIST.index("N")], "N"]
                    )
                    shift_2.append(
                        [sequence[i - 1], i, chem_shifts[i - 1][ATOM_LIST.index("C")], "C"]
                    )
                    shift_3.append(
                        [sequence[i - 1], i, chem_shifts[i - 1][c_ind], ATOM_LIST[c_ind]]
                    )
    if dimensionality == 2:
        write_2d(shift_1, shift_3, "ncocx_2d")
    elif dimensionality == 3:
        write_3d(shift_1, shift_2, shift_3, "ncocx_3d")


def ncoca_cb(sequence, chem_shifts, dimensionality):
    # this function produces peak list for inter residue NCOCX experiment
    shift_1 = []
    shift_2 = []
    shift_3 = []
    for i, residue in enumerate(sequence):  # loop through the sequence
        if (
            (i != 0)
            and (chem_shifts[i][ATOM_LIST.index("N")] != 0)
            and (chem_shifts[i - 1][ATOM_LIST.index("C")] != 0)
        ):  # check if the chemical shifts of the two are not zero
            for c_ind in c_atm_ind:
                if (chem_shifts[i - 1][c_ind] != 0) and (
                    ATOM_LIST[c_ind] == "CA" or ATOM_LIST[c_ind] == "CB"
                ):
                    shift_1.append(
                        [residue, i + 1, chem_shifts[i][ATOM_LIST.index("N")], "N"]
                    )
                    shift_2.append(
                        [sequence[i - 1], i, chem_shifts[i - 1][ATOM_LIST.index("C")], "C"]
                    )
                    shift_3.append(
                        [sequence[i - 1], i, chem_shifts[i - 1][c_ind], ATOM_LIST[c_ind]]
                    )
    if dimensionality == 2:
        write_2d(shift_1, shift_3, "ncocacb_2d")
    elif dimensionality == 3:
        write_3d(shift_1, shift_2, shift_3, "ncocacb_3d")


def dqsqsq_inter(sequence, chem_shifts):
    shift_1 = []
    shift_2 = []
    shift_3 = []
    neighbors = [
        ["C", "A"],
        ["A", "B"],
        ["B", "G"],
        ["G", "D"],
        ["D", "E"],
        ["E", "Z"],
        ["H", "Z"],
        ["E", "H"],
    ]
    for i, residue in enumerate(sequence):  # loop through the sequence
        for j, atm_pos in enumerate(simple_atom_positions):
            for k in range(6):
                if (
                    (atom_type[j] == "C")
                    and (j + 1 + k < len(atom_positions))
                    and ([atm_pos, simple_atom_positions[j + 1 + k]] in neighbors)
                ):
                    if (
                        (chem_shifts[i][j] != 0)
                        and (chem_shifts[i][j + 1 + k] != 0)
                        and (dq_exceptions(residue, ATOM_LIST[j], ATOM_LIST[j + 1 + k]) == 1)
                    ):
                        shift_1.append(
                            [
                                residue,
                                i + 1,
                                chem_shifts[i][j] + chem_shifts[i][j + 1 + k],
                                "%s+%s" % (ATOM_LIST[j], ATOM_LIST[j + 1 + k]),
                            ]
                        )
                        shift_2.append(
                            [residue, i + 1, chem_shifts[i][j], "%s" % (ATOM_LIST[j])]
                        )
                        shift_3.append(
                            [
                                residue,
                                i + 1,
                                chem_shifts[i][j + 1 + k],
                                "%s" % (ATOM_LIST[j + 1 + k]),
                            ]
                        )
                        shift_1.append(
                            [
                                residue,
                                i + 1,
                                chem_shifts[i][j + 1 + k] + chem_shifts[i][j],
                                "%s+%s" % (ATOM_LIST[j + 1 + k], ATOM_LIST[j]),
                            ]
                        )
                        shift_2.append(
                            [
                                residue,
                                i + 1,
                                chem_shifts[i][j + 1 + k],
                                "%s" % (ATOM_LIST[j + 1 + k]),
                            ]
                        )
                        shift_3.append(
                            [residue, i + 1, chem_shifts[i][j], "%s" % (ATOM_LIST[j])]
                        )
                        shift_1.append(
                            [
                                residue,
                                i + 1,
                                chem_shifts[i][j] + chem_shifts[i][j + 1 + k],
                                "%s+%s" % (ATOM_LIST[j], ATOM_LIST[j + 1 + k]),
                            ]
                        )
                        shift_2.append(
                            [residue, i + 1, chem_shifts[i][j], "%s" % (ATOM_LIST[j])]
                        )
                        shift_3.append(
                            [residue, i + 1, chem_shifts[i][j], "%s" % (ATOM_LIST[j])]
                        )
                        shift_1.append(
                            [
                                residue,
                                i + 1,
                                chem_shifts[i][j + 1 + k] + chem_shifts[i][j],
                                "%s+%s" % (ATOM_LIST[j + 1 + k], ATOM_LIST[j]),
                            ]
                        )
                        shift_2.append(
                            [
                                residue,
                                i + 1,
                                chem_shifts[i][j + 1 + k],
                                "%s" % (ATOM_LIST[j + 1 + k]),
                            ]
                        )
                        shift_3.append(
                            [
                                residue,
                                i + 1,
                                chem_shifts[i][j + 1 + k],
                                "%s" % (ATOM_LIST[j + 1 + k]),
                            ]
                        )
                        for other_f3_atm_ind in range(len(atom_positions)):
                            if (
                                (chem_shifts[i][other_f3_atm_ind] != 0)
                                and (atom_type[other_f3_atm_ind] == "C")
                                and (other_f3_atm_ind != j + 1 + k)
                                and (other_f3_atm_ind != j)
                            ):
                                shift_1.append(
                                    [
                                        residue,
                                        i + 1,
                                        chem_shifts[i][j] + chem_shifts[i][j + 1 + k],
                                        "%s+%s" % (ATOM_LIST[j], ATOM_LIST[j + 1 + k]),
                                    ]
                                )
                                shift_2.append(
                                    [
                                        residue,
                                        i + 1,
                                        chem_shifts[i][j],
                                        "%s" % (ATOM_LIST[j]),
                                    ]
                                )
                                shift_3.append(
                                    [
                                        residue,
                                        i + 1,
                                        chem_shifts[i][other_f3_atm_ind],
                                        "%s" % (ATOM_LIST[other_f3_atm_ind]),
                                    ]
                                )
                                shift_1.append(
                                    [
                                        residue,
                                        i + 1,
                                        chem_shifts[i][j + 1 + k] + chem_shifts[i][j],
                                        "%s+%s" % (ATOM_LIST[j + 1 + k], ATOM_LIST[j]),
                                    ]
                                )
                                shift_2.append(
                                    [
                                        residue,
                                        i + 1,
                                        chem_shifts[i][j + 1 + k],
                                        "%s" % (ATOM_LIST[j + 1 + k]),
                                    ]
                                )
                                shift_3.append(
                                    [
                                        residue,
                                        i + 1,
                                        chem_shifts[i][other_f3_atm_ind],
                                        "%s" % (ATOM_LIST[other_f3_atm_ind]),
                                    ]
                                )
                        if i != 0:
                            for other_f3_atm_ind in range(len(atom_positions)):
                                if (chem_shifts[i - 1][other_f3_atm_ind] != 0) and (
                                    atom_type[other_f3_atm_ind] == "C"
                                ):
                                    shift_1.append(
                                        [
                                            residue,
                                            i + 1,
                                            chem_shifts[i][j]
                                            + chem_shifts[i][j + 1 + k],
                                            "%s+%s" % (ATOM_LIST[j], ATOM_LIST[j + 1 + k]),
                                        ]
                                    )
                                    shift_2.append(
                                        [
                                            residue,
                                            i + 1,
                                            chem_shifts[i][j],
                                            "%s" % (ATOM_LIST[j]),
                                        ]
                                    )
                                    shift_3.append(
                                        [
                                            sequence[i - 1],
                                            i,
                                            chem_shifts[i - 1][other_f3_atm_ind],
                                            "%s" % (ATOM_LIST[other_f3_atm_ind]),
                                        ]
                                    )
                                    shift_1.append(
                                        [
                                            residue,
                                            i + 1,
                                            chem_shifts[i][j + 1 + k]
                                            + chem_shifts[i][j],
                                            "%s+%s" % (ATOM_LIST[j + 1 + k], ATOM_LIST[j]),
                                        ]
                                    )
                                    shift_2.append(
                                        [
                                            residue,
                                            i + 1,
                                            chem_shifts[i][j + 1 + k],
                                            "%s" % (ATOM_LIST[j + 1 + k]),
                                        ]
                                    )
                                    shift_3.append(
                                        [
                                            sequence[i - 1],
                                            i,
                                            chem_shifts[i - 1][other_f3_atm_ind],
                                            "%s" % (ATOM_LIST[other_f3_atm_ind]),
                                        ]
                                    )
                        if i != len(sequence) - 1:
                            for other_f3_atm_ind in range(len(atom_positions)):
                                if (chem_shifts[i + 1][other_f3_atm_ind] != 0) and (
                                    atom_type[other_f3_atm_ind] == "C"
                                ):
                                    shift_1.append(
                                        [
                                            residue,
                                            i + 1,
                                            chem_shifts[i][j]
                                            + chem_shifts[i][j + 1 + k],
                                            "%s+%s" % (ATOM_LIST[j], ATOM_LIST[j + 1 + k]),
                                        ]
                                    )
                                    shift_2.append(
                                        [
                                            residue,
                                            i + 1,
                                            chem_shifts[i][j],
                                            "%s" % (ATOM_LIST[j]),
                                        ]
                                    )
                                    shift_3.append(
                                        [
                                            sequence[i + 1],
                                            i + 2,
                                            chem_shifts[i + 1][other_f3_atm_ind],
                                            "%s" % (ATOM_LIST[other_f3_atm_ind]),
                                        ]
                                    )
                                    shift_1.append(
                                        [
                                            residue,
                                            i + 1,
                                            chem_shifts[i][j + 1 + k]
                                            + chem_shifts[i][j],
                                            "%s+%s" % (ATOM_LIST[j + 1 + k], ATOM_LIST[j]),
                                        ]
                                    )
                                    shift_2.append(
                                        [
                                            residue,
                                            i + 1,
                                            chem_shifts[i][j + 1 + k],
                                            "%s" % (ATOM_LIST[j + 1 + k]),
                                        ]
                                    )
                                    shift_3.append(
                                        [
                                            sequence[i + 1],
                                            i + 2,
                                            chem_shifts[i + 1][other_f3_atm_ind],
                                            "%s" % (ATOM_LIST[other_f3_atm_ind]),
                                        ]
                                    )
    write_3d(shift_1, shift_2, shift_3, "dqsqsq_inter")


def dqsqsq_intra(sequence, chem_shifts):
    shift_1 = []
    shift_2 = []
    shift_3 = []
    neighbors = [
        ["C", "A"],
        ["A", "B"],
        ["B", "G"],
        ["G", "D"],
        ["D", "E"],
        ["E", "Z"],
        ["H", "Z"],
        ["E", "H"],
    ]
    for i, residue in enumerate(sequence):  # loop through the sequence
        for j, atm_pos in enumerate(simple_atom_positions):
            for k in range(6):
                if (
                    (atom_type[j] == "C")
                    and (j + 1 + k < len(atom_positions))
                    and ([atm_pos, simple_atom_positions[j + 1 + k]] in neighbors)
                ):
                    if (chem_shifts[i][j] != 0) and (chem_shifts[i][j + 1 + k] != 0):
                        shift_1.append(
                            [
                                residue,
                                i + 1,
                                chem_shifts[i][j] + chem_shifts[i][j + 1 + k],
                                "%s+%s" % (ATOM_LIST[j], ATOM_LIST[j + 1 + k]),
                            ]
                        )
                        shift_2.append(
                            [residue, i + 1, chem_shifts[i][j], "%s" % (ATOM_LIST[j])]
                        )
                        shift_3.append(
                            [
                                residue,
                                i + 1,
                                chem_shifts[i][j + 1 + k],
                                "%s" % (ATOM_LIST[j + 1 + k]),
                            ]
                        )
                        shift_1.append(
                            [
                                residue,
                                i + 1,
                                chem_shifts[i][j + 1 + k] + chem_shifts[i][j],
                                "%s+%s" % (ATOM_LIST[j + 1 + k], ATOM_LIST[j]),
                            ]
                        )
                        shift_2.append(
                            [
                                residue,
                                i + 1,
                                chem_shifts[i][j + 1 + k],
                                "%s" % (ATOM_LIST[j + 1 + k]),
                            ]
                        )
                        shift_3.append(
                            [residue, i + 1, chem_shifts[i][j], "%s" % (ATOM_LIST[j])]
                        )
    write_3d(shift_1, shift_2, shift_3, "dqsqsq_intra")


def dqsq(sequence, chem_shifts):
    shift_1 = []
    shift_2 = []
    neighbors = [
        ["C", "A"],
        ["A", "B"],
        ["B", "G"],
        ["G", "D"],
        ["D", "E"],
        ["E", "Z"],
        ["H", "Z"],
        ["E", "H"],
    ]
    for i, residue in enumerate(sequence):  # loop through the sequence
        for j, atm_pos in enumerate(simple_atom_positions):
            for k in range(6):
                if (
                    (atom_type[j] == "C")
                    and (j + 1 + k < len(atom_positions))
                    and ([atm_pos, simple_atom_positions[j + 1 + k]] in neighbors)
                ):
                    if (chem_shifts[i][j] != 0) and (chem_shifts[i][j + 1 + k] != 0):
                        shift_1.append(
                            [
                                residue,
                                i + 1,
                                chem_shifts[i][j] + chem_shifts[i][j + 1 + k],
                                "%s+%s" % (ATOM_LIST[j], ATOM_LIST[j + 1 + k]),
                            ]
                        )
                        shift_2.append(
                            [residue, i + 1, chem_shifts[i][j], "%s" % (ATOM_LIST[j])]
                        )
                        shift_1.append(
                            [
                                residue,
                                i + 1,
                                chem_shifts[i][j + 1 + k] + chem_shifts[i][j],
                                "%s+%s" % (ATOM_LIST[j + 1 + k], ATOM_LIST[j]),
                            ]
                        )
                        shift_2.append(
                            [
                                residue,
                                i + 1,
                                chem_shifts[i][j + 1 + k],
                                "%s" % (ATOM_LIST[j + 1 + k]),
                            ]
                        )
    write_2d(shift_1, shift_2, "dqsq")


def canh(sequence, chem_shifts, dimensionality):
    # this function produces the peak list for proton detected CA-N-H experiment
    shift_1 = []
    shift_2 = []
    shift_3 = []
    for i, residue in enumerate(sequence):  # loop through the sequence
        if (
            (chem_shifts[i][ATOM_LIST.index("H")] != 0)
            and (chem_shifts[i][ATOM_LIST.index("N")] != 0)
            and (chem_shifts[i][ATOM_LIST.index("CA")] != 0)
        ):  # check if the chemical shifts of the CA, HA and N are not zero
            shift_1.append([residue, i + 1, chem_shifts[i][ATOM_LIST.index("CA")], "CA"])
            shift_2.append([residue, i + 1, chem_shifts[i][ATOM_LIST.index("N")], "N"])
            shift_3.append([residue, i + 1, chem_shifts[i][ATOM_LIST.index("H")], "H"])
    if dimensionality == 2:
        write_2d(shift_1, shift_3, "canh_2d")
    elif dimensionality == 3:
        write_3d(shift_1, shift_2, shift_3, "canh_3d")


def conh(sequence, chem_shifts, dimensionality):
    # this function produces the peak list for proton detected C-N(n+1)-HA(n+1)
    #  experiment
    shift_1 = []
    shift_2 = []
    shift_3 = []
    for i, residue in enumerate(sequence):  # loop through the sequence
        if i != 0:
            if (
                (chem_shifts[i - 1][ATOM_LIST.index("H")] != 0)
                and (chem_shifts[i - 1][ATOM_LIST.index("N")] != 0)
                and (chem_shifts[i][ATOM_LIST.index("C")] != 0)
            ):  # check if the chemical shifts of the Co, H and N are not zero
                shift_1.append([residue, i + 1, chem_shifts[i][ATOM_LIST.index("C")], "C"])
                shift_2.append(
                    [sequence[i - 1], i, chem_shifts[i - 1][ATOM_LIST.index("N")], "N"]
                )
                shift_3.append(
                    [sequence[i - 1], i, chem_shifts[i - 1][ATOM_LIST.index("H")], "H"]
                )
    if dimensionality == 2:
        write_2d(shift_1, shift_3, "conh_2d")
    elif dimensionality == 3:
        write_3d(shift_1, shift_2, shift_3, "conh_3d")


def caconh(sequence, chem_shifts, dimensionality):
    # this function produces the peak list for proton detected CA-(CO)-N-H experiment
    shift_1 = []
    shift_2 = []
    shift_3 = []
    for i, residue in enumerate(sequence):  # loop through the sequence
        if i != 0:
            if (
                (chem_shifts[i][ATOM_LIST.index("H")] != 0)
                and (chem_shifts[i][ATOM_LIST.index("N")] != 0)
                and (chem_shifts[i - 1][ATOM_LIST.index("CA")] != 0)
                and (chem_shifts[i - 1][ATOM_LIST.index("C")] != 0)
            ):
                shift_1.append(
                    [sequence[i - 1], i, chem_shifts[i - 1][ATOM_LIST.index("CA")], "CA"]
                )
                shift_2.append([residue, i + 1, chem_shifts[i][ATOM_LIST.index("N")], "N"])
                shift_3.append([residue, i + 1, chem_shifts[i][ATOM_LIST.index("H")], "H"])
    if dimensionality == 2:
        write_2d(shift_1, shift_3, "caconh")
    elif dimensionality == 3:
        write_3d(shift_1, shift_2, shift_3, "caconh")


def cocanh(sequence, chem_shifts, dimensionality):
    # this function produces the peak list for proton detected CO-(CA)-N-H experiment
    shift_1 = []
    shift_2 = []
    shift_3 = []
    for i, residue in enumerate(sequence):  # loop through the sequence
        if i != 0:
            if (
                (chem_shifts[i - 1][ATOM_LIST.index("H")] != 0)
                and (chem_shifts[i - 1][ATOM_LIST.index("N")] != 0)
                and (chem_shifts[i][ATOM_LIST.index("C")] != 0)
                and (chem_shifts[i][ATOM_LIST.index("CA")] != 0)
            ):
                shift_1.append([residue, i + 1, chem_shifts[i][42], "C"])
                shift_2.append([sequence[i - 1], i, chem_shifts[i - 1][33], "N"])
                shift_3.append([sequence[i - 1], i, chem_shifts[i - 1][0], "H"])
    if dimensionality == 2:
        write_2d(shift_1, shift_3, "cocanh_2d")
    elif dimensionality == 3:
        write_3d(shift_1, shift_2, shift_3, "cocanh_3d")


def ncah(sequence, chem_shifts, dimensionality):
    # this function produces the peak list for proton detected N-CA-H experiment
    shift_1 = []
    shift_2 = []
    shift_3 = []
    for i, residue in enumerate(sequence):  # loop through the sequence
        if i != 0:
            if (
                (chem_shifts[i][ATOM_LIST.index("H")] != 0)
                and (chem_shifts[i][ATOM_LIST.index("N")] != 0)
                and (chem_shifts[i][ATOM_LIST.index("HA")] != 0)
            ):
                shift_1.append([residue, i + 1, chem_shifts[i][ATOM_LIST.index("N")], "N"])
                shift_2.append(
                    [residue, i + 1, chem_shifts[i][ATOM_LIST.index("CA")], "CA"]
                )
                shift_3.append(
                    [residue, i + 1, chem_shifts[i][ATOM_LIST.index("HA")], "HA"]
                )
    if dimensionality == 2:
        write_2d(shift_1, shift_3, "ncah_2d")
    elif dimensionality == 3:
        write_3d(shift_1, shift_2, shift_3, "ncah_3d")


def sqsqsq_inter(sequence, chem_shifts):
    # this function produces peak list for intra-residue CC spin
    #  diffusion experiment where
    shift_1 = []
    shift_2 = []
    shift_3 = []
    for i, residue in enumerate(sequence):  # loop through the sequence
        for j, atom in enumerate(atom_type):  # loop through all atomtypes
            if atom == "C":  # check if it is a carbon atom
                for c_ind in c_atm_ind:  # loop through all the carbon indices
                    if (
                        c_ind != j
                    ):  # ensure that the two indices "j" and "c_ind" are not the same
                        for (
                            c_ind_1
                        ) in c_atm_ind:  # loop through all the carbon indices for i-1
                            if i != 0:
                                if (
                                    chem_shifts[i - 1][c_ind_1] != 0
                                    and chem_shifts[i][c_ind] != 0
                                    and chem_shifts[i][j] != 0
                                ):  # check if the chemical shifts of the three are
                                    #  not zero
                                    shift_1.append(
                                        [residue, i + 1, chem_shifts[i][j], ATOM_LIST[j]]
                                    )
                                    shift_2.append(
                                        [
                                            residue,
                                            i + 1,
                                            chem_shifts[i][c_ind],
                                            ATOM_LIST[c_ind],
                                        ]
                                    )
                                    shift_3.append(
                                        [
                                            sequence[i - 1],
                                            i,
                                            chem_shifts[i - 1][c_ind_1],
                                            ATOM_LIST[c_ind_1],
                                        ]
                                    )
                        for (
                            c_ind_1
                        ) in c_atm_ind:  # loop through all the carbon indices for i+1
                            if i != len(sequence) - 1:
                                if (
                                    chem_shifts[i + 1][c_ind_1] != 0
                                    and chem_shifts[i][c_ind] != 0
                                    and chem_shifts[i][j] != 0
                                ):  # check if the CS are not zero
                                    shift_1.append(
                                        [residue, i + 1, chem_shifts[i][j], ATOM_LIST[j]]
                                    )
                                    shift_2.append(
                                        [
                                            residue,
                                            i + 1,
                                            chem_shifts[i][c_ind],
                                            ATOM_LIST[c_ind],
                                        ]
                                    )
                                    shift_3.append(
                                        [
                                            sequence[i + 1],
                                            i + 2,
                                            chem_shifts[i + 1][c_ind_1],
                                            ATOM_LIST[c_ind_1],
                                        ]
                                    )
    write_3d(shift_1, shift_2, shift_3, "sqsqsq_inter")


def peaks_proton_heavy(sequence, chem_shifts, atm_1, atm_2, direct):
    # this function produces peak list for CH, HC, HN and NH type experiments
    # atm_2 is by default the heavy atom
    # direct defines the atom on the direct dimension
    shift_1 = []
    shift_2 = []
    for i, residue in enumerate(sequence):  # loop through the sequence
        for j, atom in enumerate(atom_type):  # loop through all atomtypes
            if atom == atm_2:  # check if it is a heavy atom
                for h_ind in h_atm_ind:  # loop on proton indices
                    if (
                        atom_positions[h_ind] == atom_positions[j]
                        and chem_shifts[i][h_ind] != 0
                        and chem_shifts[i][j] != 0
                    ):
                        shift_1.append(
                            [residue, i + 1, chem_shifts[i][h_ind], ATOM_LIST[h_ind]]
                        )
                        shift_2.append([residue, i + 1, chem_shifts[i][j], ATOM_LIST[j]])
    if atm_1 == direct:
        extension = "%s%s" % (atm_1, atm_2)
        write_2d(shift_1, shift_2, extension.lower())
    else:
        extension = "%s%s" % (atm_2, atm_1)
        write_2d(shift_2, shift_1, extension.lower())


def dq_exceptions(residue, atm_1, atm_2):
    if residue == "H":
        not_neighbors = [["CD2", "CE1"]]
        if ([atm_1, atm_2] in not_neighbors) or ([atm_2, atm_1] in not_neighbors):
            return 0
    if (residue == "F") or (residue == "Y"):
        not_neighbors = [["CD1", "CE2"], ["CD2", "CE1"]]
        if ([atm_1, atm_2] in not_neighbors) or ([atm_2, atm_1] in not_neighbors):
            return 0
    if residue == "I":
        not_neighbors = [["CD1", "CG2"]]
        if ([atm_1, atm_2] in not_neighbors) or ([atm_2, atm_1] in not_neighbors):
            return 0
    if residue == "W":
        not_neighbors = [["CD1", "CE2"], ["CZ3", "CE2"], ["CZ2", "CE3"]]
        if ([atm_1, atm_2] in not_neighbors) or ([atm_2, atm_1] in not_neighbors):
            return 0
    else:
        return 1
