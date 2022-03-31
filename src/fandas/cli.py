# NMR Spectroscopy Group, Utrecht University
import argparse
import logging

# import os
import sys

# from pathlib import Path

# import numpy as np
# from fandas.modules import experiment

from fandas.modules.experiment import Experiment

from fandas.modules.chemical import (
    STANDARD_DATA,
    #     amino_acids,
    #     list_2d,
    #     list_2dd,
    #     list_3d,
    #     list_3dd,
    #     secondary_structures,
)
from fandas.modules.experiments import (
    caconh,
    canh,
    cc_spindiff_dist,
    cc_spindiff_inter,
    cc_spindiff_intra,
    chh_dist,
    chh_dist_3d,
    chhc_dist,
    cocanh,
    conh,
    dqsq,
    dqsqsq_inter,
    dqsqsq_intra,
    hh,
    hh_dist,
    nca,
    ncacx,
    ncacx_dist,
    ncacx_inter,
    ncah,
    nco,
    ncoca_cb,
    ncocx,
    ncocx_dist,
    nhh_dist,
    nhh_dist_3d,
    nhhc_dist,
    peaks_proton_heavy,
    sqsqsq_inter,
)
from fandas.modules.input import InputFile
from fandas.modules.utils import (
    dict2array,
    assign_chemical_shifts,
    check_user_input,
    fractional_deuteration,
    glycerol_label,
    replace_bmrb,
    rev_label,
    load_sequence,
    load_secondary_structure,
    glycerol_label_13,
)
from fandas.modules.chemical_shift import ChemShift
from fandas.version import VERSION

# Setup logging
log = logging.getLogger("fandaslog")
ch = logging.StreamHandler()
formatter = logging.Formatter(
    " [%(asctime)s %(funcName)s:L%(lineno)d %(levelname)s] %(message)s"
)
ch.setFormatter(formatter)
log.addHandler(ch)

__author__ = ["Siddarth Narasimhan", "Rodrigo Honorato"]

# ===========================================================================================================
# Define arguments

ap = argparse.ArgumentParser()
ap.add_argument(
    "input_file",
    help="",
)
# ap.add_argument(
#     "-seq",
#     help="input sequence as a text file (REQUIRED)",
#     required=True,
#     dest="sequence_fname",
# )
# ap.add_argument(
#     "-ss",
#     help=(
#         "secondary structures as a text file: 'a'- alpha helix, 'b'- beta sheet &"
#         " 'c'- random coil; if unspecified, will use BMRB averages: 'n'"
#     ),
#     dest="secondary_structure_fname",
# )

# ap.add_argument(
#     "-bt",
#     help="BMRB tables as a text file in NMR Star format",
#     default="",
#     dest="bmrb_table_fname",
# )

# ap.add_argument(
#     "-bt_seq_start",
#     help="Which should be the first residue number of your sequence.",
#     type=int,
#     default=1,
# )
# ap.add_argument(
#     "-btc",
#     help=(
#         "column numbers for residue number, atom name and chemical shifts in the"
#         " BMRB tables"
#     ),
#     nargs=3,
#     default=[1, 3, 4],
# )
# ap.add_argument(
#     "-ls",
#     help=(
#         "labelling scheme. Default is Fully Labelled. Other options: fw = Forward"
#         " Labelled, rv = Reverse Labelled, gl13 = 1,3 Glycerol Labelling, gl2= 2"
#         " Glycerol Labelling"
#     ),
#     dest="labeling_scheme",
#     default="fully",
#     choices=["fully", "forward", "reverse", "gl13", "gl2"],
# )

# ap.add_argument(
#     "-dl",
#     help=(
#         "list of 13C & 15N (for forward labelling) or 12C & 14N (for reverse"
#         " labelling) amino acids"
#     ),
#     nargs="+",
#     default="",
# )
# ap.add_argument(
#     "-cl",
#     help=(
#         "list of 13C (for forward labelling) or 12C (for reverse labelling) amino"
#         " acids"
#     ),
#     nargs="+",
#     default="",
# )
# ap.add_argument(
#     "-nl",
#     help=(
#         "list of 15N (for forward labelling) or 14N (for reverse labelling) amino"
#         " acids"
#     ),
#     nargs="+",
#     default="",
# )
# ap.add_argument(
#     "-fd",
#     help=(
#         "fractionally deuterated, if you use this flag, it would be implemented"
#         " automatically"
#     ),
#     action="store_true",
#     dest="fractionally_deuterated",
# )
# ap.add_argument(
#     "-exp_2d", help=f"List of 2D experiments: {list_2d}", nargs="+", default=[]
# )
# ap.add_argument(
#     "-exp_2dd",
#     help=(
#         "List of distance encoded (distance list, -dl and limit -dlm  must be"
#         f" provided) 2D experiments: {list_2dd}"
#     ),
#     nargs="+",
#     default=[],
# )
# ap.add_argument(
#     "-exp_3d", help=f"List of 3D experiments: {list_3d}", nargs="+", default=[]
# )

# ap.add_argument(
#     "-exp_3dd",
#     help=(
#         "List of distance encoded (distance list, -dl and limit -dlm must be"
#         f" provided) 3D experiments: {list_3dd}"
#     ),
#     nargs="+",
#     default=[],
# )
# ap.add_argument(
#     "-sl",
#     help="automatically assign peak labels in the sparky file",
#     action="store_true",
# )

# ap.add_argument(
#     "-dlist",
#     help=(
#         "Distance list file in angstorms containing "
#         "resi_num,atm_nam,resi_num_2,atm_nam_2,dist - EXAMPLE: 2,CA,4,CB,4.5"
#     ),
#     default="",
#     dest="distance_list_fname",
# )
# ap.add_argument(
#     "-dlim",
#     help="Distance limit in angstorms",
#     default=5,
#     dest="distance_limit",
# )
# ap.add_argument(
#     "-o",
#     help="names for output",
#     default="fandas_output",
#     dest="output_fname",
# )

ap.add_argument(
    "-v",
    "--version",
    help="show version",
    action="version",
    version=f"Running {ap.prog} v{VERSION}",
)


def load_args(ap):
    """Load argument parser"""
    return ap.parse_args()


# ====================================================================================#
# Define CLI


def cli(ap, main):
    """Command-line interface entry point."""
    cmd = load_args(ap)
    main(**vars(cmd))


def maincli():
    """Execute main client."""
    cli(ap, main)


# ====================================================================================#
# Main code


def main(
    input_file,
    # bmrb_table_fname,
    # bt_seq_start,
    # btc,
    # cl,
    # dl,
    # distance_limit,
    # distance_list_fname,
    # exp_2d,
    # exp_2dd,
    # exp_3d,
    # exp_3dd,
    # fractionally_deuterated,
    # labeling_scheme,
    # nl,
    # output_fname,
    # sequence_fname,
    # sl,
    # secondary_structure_fname,
    log_level="DEBUG",
):

    # Start #=========================================================================#
    log.setLevel(log_level)
    log.info("########################################")
    log.info(f"     Welcome to FANDAS {VERSION}")
    log.info("########################################")

    inp = InputFile(input_file)

    # Read input sequence # ==========================================================#
    sequence = inp.data["general"]["sequence"]

    # Assign single letter secondary structure # =====================================#

    secondary_structure = inp.data["general"]["secondary_structure"]
    if not secondary_structure:
        log.warning('No secondary structure input given, assuming "n"')
        secondary_structure = ["n"] * len(sequence)

    # Assign bmrb average to make up for missing ss. If the ss assi
    # FIXME: this does not account for GAPS
    if len(secondary_structure) < len(sequence):
        difference = len(sequence) - len(secondary_structure)
        for i in range(difference):
            secondary_structure += "n"

    # #If the user defined secondary structure is longer than sequence, chop it##
    elif len(secondary_structure) > len(sequence):
        secondary_structure = secondary_structure[: len(sequence)]

    # Assign chemical shifts #========================================================#
    log.info("Assigning average shifts")

    chem_shifts = ChemShift(sequence, secondary_structure)
    chem_shifts.assign(STANDARD_DATA)

    # Replace the average shifts with provided shifts in the form of BMRB table ======#
    bmrb_table_fname = inp.data["BMRB"]["bmrb_table_fname"]
    if bmrb_table_fname:

        resnum_col = inp.data["BMRB"]["resnum_column"]
        atom_col = inp.data["BMRB"]["atom_column"]
        shift_col = inp.data["BMRB"]["chemical_shift_column"]
        bt_seq_start = inp.data["BMRB"]["sequence_start"]

        log.info(
            "Replacing the average shifts with provided "
            "shifts in the form of BMRB table"
        )
        chem_shifts.replace_with_bmrb(
            bmrb_table_fname, resnum_col, atom_col, shift_col, bt_seq_start
        )

    # Incorporate forward, reverse & glycerol labeling schemes #=====================#
    #
    #  Reverse = remove from list (set to 0.0)
    #  Glycerol = remove based on a pre-set list
    #  Fully = Use everything
    labeling_params = inp.data["labeling"]

    # configure labelled amino acids list if fwd or rv labeling schemes are used
    log.info(f"Incorporating labeling scheme: {labeling_params['scheme']}")

    chem_shifts.label(labeling_params)

    # Take fractional deuteration in account #========================================#
    if inp.data["general"]["fractionally_deuterated"]:
        chem_shifts.consider_deuteration()

    # Read distance lists #===========================================================#

    # # ##(8) Read distance lists
    # if distance_list_fname != "":
    #     distance_list_fname = Path(distance_list_fname).resolve()
    #     log.info(f"Reading distance lists from {distance_list_fname}")
    #     distance_limit = distance_limit
    #     distances = []
    #     dist_list = []
    #     with open(distance_list_fname, "r") as list_file:
    #         for line in list_file:
    #             line = line.rstrip()
    #             line = line.replace(" ", "")
    #             line = line.split(",")
    #             if float(line[4]) <= float(distance_limit):
    #                 distances.append(line[4])
    #                 dist_list.append([int(line[0]), line[1], int(line[2]), line[3]])
    # else:
    #     dist_list = []

    # # ##(9) Round off the final chemical shifts to 2 decimal places and include
    # # options for sparky labels
    # log.info("Round off the final chemical shifts to 2 decimal places")
    # chem_shifts = np.around(chem_shifts, decimals=2)
    # global sl
    # sl = sl

    # Make predictions #==============================================================#
    experiment = Experiment(chem_shifts, inp.data)
    experiment.run_experiments()

    # if exp_3d != []:
    #     for experiment in exp_3d:
    #         experiment = experiment.upper()
    #         if experiment == "NCACX":
    #             ncacx(sequence, chem_shifts, 3)
    #         elif experiment == "NCACX_INTER":
    #             ncacx_inter(sequence, chem_shifts, 3)
    #         elif experiment == "NCOCX":
    #             ncocx(sequence, chem_shifts, 3)
    #         elif experiment == "NCOCA_CB":
    #             ncoca_cb(sequence, chem_shifts, 3)
    #         elif experiment == "SQSQSQ_INTER":
    #             sqsqsq_inter(sequence, chem_shifts)
    #         elif experiment == "DQSQSQ_INTRA":
    #             dqsqsq_intra(sequence, chem_shifts)
    #         elif experiment == "DQSQSQ_INTER":
    #             dqsqsq_inter(sequence, chem_shifts)
    #         elif experiment == "CANH":
    #             canh(sequence, chem_shifts, 3)
    #         elif experiment == "CONH":
    #             conh(sequence, chem_shifts, 3)
    #         elif experiment == "CACONH":
    #             caconh(sequence, chem_shifts, 3)
    #         elif experiment == "COCANH":
    #             conh(sequence, chem_shifts, 3)
    #         elif experiment == "NCAH":
    #             ncah(sequence, chem_shifts, 3)
    #         else:
    #             logging.warning("%s is not a valid experiment" % experiment)

    # if (exp_2dd != []) and (dist_list != []):
    #     for experiment in exp_2dd:
    #         experiment = experiment.upper()
    #         if experiment == "CC_SPINDIFF":
    #             cc_spindiff_dist(sequence, chem_shifts, dist_list, distances)
    #         elif experiment == "HH":
    #             hh_dist(sequence, chem_shifts, dist_list, distances)
    #         elif experiment == "CHHC":
    #             chhc_dist(sequence, chem_shifts, dist_list, distances)
    #         elif experiment == "NHHC":
    #             nhhc_dist(sequence, chem_shifts, dist_list, distances)
    #         elif experiment == "CHH":
    #             chh_dist(sequence, chem_shifts, dist_list, distances)
    #         elif experiment == "NHH":
    #             nhh_dist(sequence, chem_shifts, dist_list, distances)
    #         # FIXME: this is not defined
    #         # elif experiment == "HHC":
    #         #     hhc_dist(chem_shifts, dist_list, distances)
    #         elif experiment == "NCACX":
    #             ncacx_dist(sequence, chem_shifts, 2, dist_list, distances)
    #         elif experiment == "NCOCX":
    #             ncocx_dist(sequence, chem_shifts, 2, dist_list, distances)
    #         else:
    #             logging.warning("%s is not a valid experiment" % experiment)

    # if (exp_3dd != []) and (dist_list != []):
    #     for experiment in exp_3dd:
    #         if experiment == "CHH":
    #             chh_dist_3d(sequence, chem_shifts, dist_list, distances)
    #         elif experiment == "NHH":
    #             nhh_dist_3d(sequence, chem_shifts, dist_list, distances)
    #         elif experiment == "NCACX":
    #             ncacx_dist(sequence, chem_shifts, 3, dist_list, distances)
    #         elif experiment == "NCOCX":
    #             ncocx_dist(sequence, chem_shifts, 3, dist_list, distances)
    #         else:
    #             logging.warning("%s is not a valid experiment" % experiment)


if __name__ == "__main__":
    sys.exit(maincli())
