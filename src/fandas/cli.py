# NMR Spectroscopy Group, Utrecht University
import argparse
import logging
import os
import sys
from pathlib import Path

import numpy as np

from fandas.modules.chemical import (
    amino_acids,
    list_2d,
    list_2dd,
    list_3d,
    list_3dd,
    secondary_structures,
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
from fandas.modules.utils import (
    dict2array,
    assign_chemical_shifts,
    check_user_input,
    fractional_deuteration,
    glycerol_label,
    replace_bmrb,
    rev_label,
)
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
ap.add_argument("-seq", help="input sequence as a text file (REQUIRED)", required=True)
ap.add_argument(
    "-ss",
    help=(
        "secondary structures as a text file: 'a'- alpha helix, 'b'- beta sheet &"
        " 'c'- random coil; if unspecified, will use BMRB averages: 'n'"
    ),
    default="",
)
ap.add_argument("-bt", help="BMRB tables as a text file in NMR Star format", default="")
ap.add_argument(
    "-bt_seq_start",
    help="Which should be the first residue number of your sequence.",
    type=int,
    default=1,
)
ap.add_argument(
    "-btc",
    help=(
        "column numbers for residue number, atom name and chemical shifts in the"
        " BMRB tables"
    ),
    nargs=3,
    default=[1, 3, 4],
)
ap.add_argument(
    "-ls",
    help=(
        "labelling scheme. Default is Fully Labelled. Other options: fw = Forward"
        " Labelled, rv = Reverse Labelled, gl13 = 1,3 Glycerol Labelling, gl2= 2"
        " Glycerol Labelling"
    ),
)
ap.add_argument(
    "-dl",
    help=(
        "list of 13C & 15N (for forward labelling) or 12C & 14N (for reverse"
        " labelling) amino acids"
    ),
    nargs="+",
    default="",
)
ap.add_argument(
    "-cl",
    help=(
        "list of 13C (for forward labelling) or 12C (for reverse labelling) amino"
        " acids"
    ),
    nargs="+",
    default="",
)
ap.add_argument(
    "-nl",
    help=(
        "list of 15N (for forward labelling) or 14N (for reverse labelling) amino"
        " acids"
    ),
    nargs="+",
    default="",
)
ap.add_argument(
    "-fd",
    help=(
        "fractionally deuterated, if you use this flag, it would be implemented"
        " automatically"
    ),
    action="store_true",
)
ap.add_argument(
    "-exp_2d", help="list of 2D experiments: %s" % list_2d, nargs="+", default=""
)
ap.add_argument(
    "-exp_2dd",
    help=(
        "list of distance encoded (distance list, -dl and limit -dlm  must be"
        " provided) 2D experiments: %s"
    )
    % list_2dd,
    nargs="+",
    default="",
)
ap.add_argument(
    "-exp_3d", help="list of 3D experiments: %s" % list_3d, nargs="+", default=""
)
ap.add_argument(
    "-exp_3dd",
    help=(
        "list of distance encoded (distance list, -dl and limit -dlm must be"
        " provided) 3D experiments: %s"
    )
    % list_3dd,
    nargs="+",
    default="",
)
ap.add_argument(
    "-sl",
    help="automatically assign peak labels in the sparky file",
    action="store_true",
)
ap.add_argument(
    "-dlist",
    help=(
        "distance list (as a file) in angstorms as follows:"
        ' "resi_num,atm_nam,resi_num_2,atm_nam_2,dist". EXAMPLE:"2,CA,4,CB,4.5"'
    ),
    default="",
)
ap.add_argument(
    "-dlim", help="distance limit in angstorms, default: 5 Angstorm", default=5
)
ap.add_argument(
    "-o", help='names for output, default: "fandas_output"', default="fandas_output"
)

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


# ===========================================================================================================
# Define CLI


def cli(ap, main):
    """Command-line interface entry point."""
    cmd = load_args(ap)
    main(**vars(cmd))


def maincli():
    """Execute main client."""
    cli(ap, main)


# ===========================================================================================================
# Main code


def main(
    bt,
    bt_seq_start,
    btc,
    cl,
    dl,
    dlim,
    dlist,
    exp_2d,
    exp_2dd,
    exp_3d,
    exp_3dd,
    fd,
    ls,
    nl,
    o,
    seq,
    sl,
    ss,
    log_level="INFO",
):
    log.setLevel(log_level)
    log.info("########################################")
    log.info(f"     Welcome to FANDAS {VERSION}")
    log.info("########################################")
    # ##(1) give sensible names to form inputs and parse inputs
    global project_name
    # global sequence
    # ##STEP 1: READ INPUT SEQUENCE AND CHECK THE
    #  SEQUENCE USING check_user_input FUNCTION# ##
    log.info(f"Reading sequence from {seq}")
    project_name = o
    sequence = Path(seq).resolve()
    if not sequence.exists():
        log.error(f"{seq} not found.")
        sys.exit()
    with open(sequence, "r") as input_fh:
        sequence = "".join(input_fh.readlines())

    # remove \n
    sequence = sequence.split(os.linesep)[0]
    # remove spaces
    sequence = "".join(sequence.split())
    # make it uppercase
    sequence = sequence.upper()
    check_user_input(sequence, amino_acids, "Sequence")

    # ##STEP 2: ASSIGN SINGLE LETTER SECONDARY STRUCTURE FOR THE RESIDUES
    # #If the secondary structure is not given, BMRB average shifts
    #   (as on 21/07/16) will be assigned##
    if ss == "":
        log.warning('No secondary structure input given, assuming "n"')
        sec_struc = ["n"] * len(sequence)
    else:
        log.info(f"Assigning secondary structure {ss}")
        ss = Path(ss).resolve()
        if not ss.exists():
            log.error(f"Secondary structure file {ss} does not exist")
            sys.exit()
        with open(ss, "r") as ss_file:
            sec_struc = list(ss_file.readline())
        for ss in sec_struc:
            if ss == "\n":
                sec_struc.remove("\n")
            if ss == " ":
                sec_struc.remove(" ")
        for i in range(len(sec_struc)):
            sec_struc[i] = sec_struc[i].lower()
        # #Assign bmrb average to make up for missing ss. If the ss assi
        if len(sec_struc) < len(sequence):
            difference = len(sequence) - len(sec_struc)
            for i in range(difference):
                sec_struc.append("n")
        # #If the user defined secondary structure is longer than sequence, chop it##
        elif len(sec_struc) > len(sequence):
            sec_struc = sec_struc[: len(sequence)]
        check_user_input(sec_struc, secondary_structures, "Secondary Structure")

    # ##(4) assign chemical shifts
    log.info("Assigning average shifts")
    chem_shifts = assign_chemical_shifts(sequence, sec_struc)

    # ##(5) replace the average shifts with provided shifts in the form of BMRB table
    if bt != "":
        log.info(
            "Replacing the average shifts with provided "
            "shifts in the form of BMRB table"
        )
        chem_shifts = replace_bmrb(sequence, chem_shifts, bt, btc, bt_seq_start)

    # convert chemical shift to a numpy array
    chem_shifts = dict2array(chem_shifts)

    # ##(6) incorporate forward, reverse & glycerol labelling schemes
    # configure labelled amino acids list if fwd or rv labelling schemes are used
    log.info(f"Incorporating labelling scheme: {ls}")
    labelling_scheme = ls
    dl = dl
    cl = cl
    nl = nl
    if (labelling_scheme == "fw") or (labelling_scheme == "rv"):
        fdl = []
        fcl = []
        fnl = []
        rdl = []
        rcl = []
        rnl = []
        if dl != []:
            if labelling_scheme == "fw":
                fdl = dl
            else:
                rdl = dl
        if cl != []:
            if labelling_scheme == "fw":
                fcl = cl
            else:
                rcl = cl
        if nl != []:
            if labelling_scheme == "fw":
                fnl = nl
            else:
                rnl = nl
        # clone amino_acid list onto aa_rm from which all the fwd_labelled amino
        #  acids will be removed
        if labelling_scheme == "fw":
            aa_rm = amino_acids
        # remove all the fwd_cn labelled amino acids from the aa_rm list
        if fdl != []:
            for fw_double in fdl:
                aa_rm.remove(fw_double.upper())
        # remove all the 13C shifts for fwd_n labelled
        # remove all the 15N labelled amino acids from the aa_rm list
        if fnl != []:
            for aa in fnl:
                chem_shifts = rev_label(sequence, chem_shifts, "rev_c", aa.upper())
                aa_rm.remove(aa.upper())
        # remove all the 15N shifts for fwd_c labelled
        # remove all the 13C labelled amino acids from the aa_rm list
        if fcl != []:
            for aa in fcl:
                chem_shifts = rev_label(sequence, chem_shifts, "rev_n", aa.upper())
                aa_rm.remove(aa.upper())
        # remove all the 13C & 15N shifts for rev_cn labelled
        if rdl != []:
            for aa in rdl:
                chem_shifts = rev_label(sequence, chem_shifts, "rev_cn", aa.upper())
        # remove all the 15N shifts for rv_n laabelled
        if rnl != []:
            for aa in rnl:
                chem_shifts = rev_label(sequence, chem_shifts, "rev_n", aa.upper())
        # remove all the 13C shifts for fwd_c labelled
        if rcl != []:
            for aa in rcl:
                chem_shifts = rev_label(sequence, chem_shifts, "rev_c", aa.upper())
        # remove all the 13C & 15N shifts for everything in the aa_rm list
        if (fdl != []) or (fcl != []) or (fnl != []):
            for aa in aa_rm:
                chem_shifts = rev_label(sequence, chem_shifts, "rev_cn", aa.upper())
    elif labelling_scheme == "gl13":
        chem_shifts = glycerol_label(sequence, chem_shifts, 13)
    elif labelling_scheme == "gl2":
        chem_shifts = glycerol_label(sequence, chem_shifts, 2)

    # ##7 TAKE FRACTIONAL DEUTERATION INTO ACCOUNT# ##
    if fd is True:
        log.info("Using fractional deuteration")
        chem_shifts = fractional_deuteration(sequence, chem_shifts)
    chem_shifts = np.around(chem_shifts, decimals=4)

    # ##(8) Read distance lists
    if dlist != "":
        dlist = Path(dlist).resolve()
        log.info(f"Reading distance lists from {dlist}")
        dlim = dlim
        distances = []
        dist_list = []
        with open(dlist, "r") as list_file:
            for line in list_file:
                line = line.rstrip()
                line = line.replace(" ", "")
                line = line.split(",")
                if float(line[4]) <= float(dlim):
                    distances.append(line[4])
                    dist_list.append([int(line[0]), line[1], int(line[2]), line[3]])
    else:
        dist_list = []

    # ##(9) Round off the final chemical shifts to 2 decimal places and include
    # options for sparky labels
    log.info("Round off the final chemical shifts to 2 decimal places")
    chem_shifts = np.around(chem_shifts, decimals=2)
    # global sl
    # sl = sl

    # ##(10) Make Predictions
    log.info("Making predictions...")
    if exp_2d:
        log.info(f"2D: {', '.join(exp_2d)}")
    if exp_2dd:
        log.info(f"2DD: {', '.join(exp_2dd)}")
    if exp_3d:
        log.info(f"3D: {', '.join(exp_3d)}")
    if exp_3dd:
        log.info(f"3DD: {', '.join(exp_3dd)}")
    if not any([exp_2d, exp_2dd, exp_3d, exp_3dd]):
        log.error("> No experiment was selected.")
        sys.exit()
    if exp_2d != []:
        for experiment in exp_2d:
            experiment = experiment.upper()
            if experiment == "NH":
                peaks_proton_heavy(sequence, chem_shifts, "H", "N", "H")
            elif experiment == "HN":
                peaks_proton_heavy(sequence, chem_shifts, "H", "N", "N")
            elif experiment == "CH":
                peaks_proton_heavy(sequence, chem_shifts, "H", "C", "C")
            elif experiment == "HC":
                peaks_proton_heavy(sequence, chem_shifts, "H", "C", "H")
            elif experiment == "HH":
                hh(sequence, chem_shifts)
            elif experiment == "DQSQ":
                dqsq(sequence, chem_shifts)
            elif experiment == "CC_SPINDIFF_INTRA":
                cc_spindiff_intra(sequence, chem_shifts)
            elif experiment == "CC_SPINDIFF_INTER":
                cc_spindiff_inter(sequence, chem_shifts)
            elif experiment == "NCA":
                nca(sequence, chem_shifts)
            elif experiment == "NCO":
                nco(sequence, chem_shifts)
            elif experiment == "NCACX":
                ncacx(sequence, chem_shifts, 2)
            elif experiment == "NCACX_INTER":
                ncacx_inter(sequence, chem_shifts, 2)
            elif experiment == "NCOCX":
                ncocx(sequence, chem_shifts, 2)
            elif experiment == "NCOCA_CB":
                ncoca_cb(sequence, chem_shifts, 2)
            elif experiment == "CANH":
                canh(sequence, chem_shifts, 2)
            elif experiment == "CONH":
                conh(sequence, chem_shifts, 2)
            elif experiment == "CACONH":
                caconh(sequence, chem_shifts, 2)
            elif experiment == "COCANH":
                cocanh(sequence, chem_shifts, 2)
            elif experiment == "NCAH":
                ncah(sequence, chem_shifts, 2)
            else:
                logging.warning("%s is not a valid experiment" % experiment)

    if exp_3d != []:
        for experiment in exp_3d:
            experiment = experiment.upper()
            if experiment == "NCACX":
                ncacx(sequence, chem_shifts, 3)
            elif experiment == "NCACX_INTER":
                ncacx_inter(sequence, chem_shifts, 3)
            elif experiment == "NCOCX":
                ncocx(sequence, chem_shifts, 3)
            elif experiment == "NCOCA_CB":
                ncoca_cb(sequence, chem_shifts, 3)
            elif experiment == "SQSQSQ_INTER":
                sqsqsq_inter(sequence, chem_shifts)
            elif experiment == "DQSQSQ_INTRA":
                dqsqsq_intra(sequence, chem_shifts)
            elif experiment == "DQSQSQ_INTER":
                dqsqsq_inter(sequence, chem_shifts)
            elif experiment == "CANH":
                canh(sequence, chem_shifts, 3)
            elif experiment == "CONH":
                conh(sequence, chem_shifts, 3)
            elif experiment == "CACONH":
                caconh(sequence, chem_shifts, 3)
            elif experiment == "COCANH":
                conh(sequence, chem_shifts, 3)
            elif experiment == "NCAH":
                ncah(sequence, chem_shifts, 3)
            else:
                logging.warning("%s is not a valid experiment" % experiment)

    if (exp_2dd != []) and (dist_list != []):
        for experiment in exp_2dd:
            experiment = experiment.upper()
            if experiment == "CC_SPINDIFF":
                cc_spindiff_dist(sequence, chem_shifts, dist_list, distances)
            elif experiment == "HH":
                hh_dist(sequence, chem_shifts, dist_list, distances)
            elif experiment == "CHHC":
                chhc_dist(sequence, chem_shifts, dist_list, distances)
            elif experiment == "NHHC":
                nhhc_dist(sequence, chem_shifts, dist_list, distances)
            elif experiment == "CHH":
                chh_dist(sequence, chem_shifts, dist_list, distances)
            elif experiment == "NHH":
                nhh_dist(sequence, chem_shifts, dist_list, distances)
            # FIXME: this is not defined
            # elif experiment == "HHC":
            #     hhc_dist(chem_shifts, dist_list, distances)
            elif experiment == "NCACX":
                ncacx_dist(sequence, chem_shifts, 2, dist_list, distances)
            elif experiment == "NCOCX":
                ncocx_dist(sequence, chem_shifts, 2, dist_list, distances)
            else:
                logging.warning("%s is not a valid experiment" % experiment)

    if (exp_3dd != []) and (dist_list != []):
        for experiment in exp_3dd:
            if experiment == "CHH":
                chh_dist_3d(sequence, chem_shifts, dist_list, distances)
            elif experiment == "NHH":
                nhh_dist_3d(sequence, chem_shifts, dist_list, distances)
            elif experiment == "NCACX":
                ncacx_dist(sequence, chem_shifts, 3, dist_list, distances)
            elif experiment == "NCOCX":
                ncocx_dist(sequence, chem_shifts, 3, dist_list, distances)
            else:
                logging.warning("%s is not a valid experiment" % experiment)


if __name__ == "__main__":
    sys.exit(maincli())
