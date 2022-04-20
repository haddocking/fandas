import argparse
import logging
import sys

from fandas.modules.chemical_shift import ChemShift
from fandas.modules.experiment import Experiment
from fandas.modules.input import InputFile
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

    # If the user defined secondary structure is longer than sequence, chop it
    elif len(secondary_structure) > len(sequence):
        secondary_structure = secondary_structure[: len(sequence)]

    # Assign chemical shifts #========================================================#
    log.info("Assigning average shifts")

    chem_shifts = ChemShift(sequence, secondary_structure)
    chem_shifts.assign()

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
    labeling_params = inp.data["labeling"]

    # configure labelled amino acids list if fwd or rv labeling schemes are used
    log.info(f"Incorporating labeling scheme: {labeling_params['scheme']}")

    chem_shifts.label(labeling_params)

    # Take fractional deuteration in account #========================================#
    if inp.data["general"]["fractionally_deuterated"]:
        chem_shifts.consider_deuteration()

    # Make predictions #==============================================================#
    experiment = Experiment(chem_shifts, inp.data)
    experiment.run_experiments()


# ====================================================================================#
if __name__ == "__main__":
    sys.exit(maincli())
