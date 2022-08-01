import logging
import sys
from functools import partial

import pandas as pd

from fandas.modules.bmrb import BMRB
from fandas.modules.chemical import (
    AA_REF,
    C_ATOMS,
    DEUTERATION,
    GL_2,
    GL_13,
    N_ATOMS,
    SS_REF,
    STANDARD_DATA,
)
from fandas.modules.residue import Residue

log = logging.getLogger("fandaslog")


class ChemShift:
    def __init__(self, sequence, secondary_structure, bmrb_table, bmrb_entity_id):
        self.residues = {}
        self.sequence = sequence
        self.secondary_structure = secondary_structure
        self.bmrb_table = bmrb_table
        self.bmrb_entity_id = bmrb_entity_id
        self.assign()

    def assign(self, data=STANDARD_DATA):
        """Assign the chemical shifts."""
        log.info(f"Loading standard chemical shifts from {data}")
        df = pd.read_csv(data)

        # atom_list = list(df.columns[2:])
        for resnum, (resname, ss) in enumerate(
            zip(self.sequence, self.secondary_structure), start=1
        ):

            # make the res one-letter
            three_letter_resname = AA_REF[resname]

            # get the full ss name
            ss = SS_REF[ss]

            sub_values = df[
                (df["RESNAME"] == three_letter_resname)
                & (df["SECONDARY_STRUCTURE"] == ss)
            ]
            atom_shift_dic = dict(
                (e, v) for e, v in zip(df.columns[2:], sub_values.iloc[0].values[2:])
            )

            self.residues[resnum] = Residue(resnum, resname, ss, atom_shift_dic)

        if self.bmrb_table:
            try:
                self.use_bmrb()
            except Exception as e:
                log.exception(e)
                sys.exit()

    def use_bmrb(self):
        """Use the BMRB table."""
        entry = BMRB(self.bmrb_table, self.bmrb_entity_id)
        log.info(f"Aligning input sequence with BMRB entity {self.bmrb_entity_id}")
        alignment_dict = entry.align_to(self.sequence)
        log.info("Updating shifts based on BMRB table")
        self._update_shifts(entry, alignment_dict)

    def _update_shifts(self, bmrb_entry, alignment_dict):
        """Update the chemical shifts based on the BMRB table."""
        for element in enumerate(self.residues.items()):
            position, _ = element
            _, standard_residue = list(self.residues.items())[position]
            matching_pos = alignment_dict[position + 1]
            _, bmrb_residue = list(bmrb_entry.residues.items())[matching_pos - 1]

            assert bmrb_residue.resname == standard_residue.resname

            for atom in standard_residue.shifts:
                if atom in bmrb_residue.shifts:
                    log.info(
                        f"# (bmrb.{bmrb_residue.resnum}) {bmrb_residue.resnum}\t"
                        f"{bmrb_residue.resname}\t{atom}\t"
                        f"{standard_residue.shifts[atom]}\t"
                        f">> {bmrb_residue.shifts[atom]}"
                    )
                    standard_residue.shifts[atom] = bmrb_residue.shifts[atom]

    def label(self, params):
        """Apply a labeling scheme."""
        labeling_func = self._get_label_func(**params)
        labeling_func()

    def _get_label_func(self, scheme, **kwargs):
        """Find out which labelling function should be used."""
        if scheme == "gl13":
            return partial(self.apply_gl13_labelling)

        elif scheme == "gl2":
            return partial(self.apply_gl2_labelling)

        elif scheme == "forward":
            return partial(
                self.apply_fw_labelling,
                fw_13c_15n=kwargs["forward"]["fw_13c_15n"],
                fw_13c=kwargs["forward"]["fw_13c"],
                fw_15n=kwargs["forward"]["fw_15n"],
            )

        elif scheme == "reverse":
            return partial(
                self.apply_rev_labelling,
                rev_12c_14n=kwargs["reverse"]["rev_12c_14n"],
                rev_12c=kwargs["reverse"]["rev_12c"],
                rev_14n=kwargs["reverse"]["rev_14n"],
            )

        elif scheme == "fully":
            return partial(self.apply_fully_labelling)

    def apply_gl13_labelling(self):
        """Label according to 1,3-13C-glycerol table (?).

        ...

        """
        self.zero_shift(GL_13)

    def apply_gl2_labelling(self):
        """Label according to 2-13C-glycerol table (?).

        ...

        """

        self.zero_shift(GL_2)

    def apply_fw_labelling(self, fw_13c_15n, fw_13c, fw_15n):
        """Apply the foward labeling scheme."""
        if fw_13c_15n:
            # remove N and C atoms
            self.zero_shift(
                C_ATOMS + N_ATOMS,
                fw_13c_15n,
            )

        if fw_13c:
            self.zero_shift(C_ATOMS, fw_13c)

        if fw_15n:
            self.zero_shift(N_ATOMS, fw_15n)

    def apply_rev_labelling(self, rev_12c_14n, rev_12c, rev_14n):
        """Apply the reverse labeling scheme."""
        # TODO: Is this the same as FW scheme??
        if rev_12c_14n:
            # remove N and C atoms
            self.zero_shift(C_ATOMS + N_ATOMS, rev_12c_14n)

        if rev_12c:
            self.zero_shift(C_ATOMS, rev_12c)

        if rev_14n:
            self.zero_shift(N_ATOMS, rev_14n)

    def apply_fully_labelling(self):
        """Apply the full labeling."""

        # Do nothing :)

        pass

    def consider_deuteration(self):
        """Consider the deuteration."""
        self.zero_shift(DEUTERATION)

    def zero_shift(self, atoms_to_be_set_to_zero, resnum_list=None):
        """Zero the shift of a residue list based on an atom iterable."""
        for resnum in self.residues:
            resname = self.residues[resnum].resname

            if resnum_list:
                if resnum in resnum_list:
                    res_check = True
                else:
                    res_check = False
            else:
                res_check = True

            if isinstance(atoms_to_be_set_to_zero, dict):
                try:
                    atom_del_list = atoms_to_be_set_to_zero[resname]
                except KeyError:
                    # no atoms to be deleted in this residue
                    continue

            elif isinstance(atoms_to_be_set_to_zero, list):
                atom_del_list = atoms_to_be_set_to_zero

            elif isinstance(atoms_to_be_set_to_zero, str):
                atom_del_list = [atoms_to_be_set_to_zero]

            if res_check:
                for atom in self.residues[resnum].shifts:
                    if atom in atom_del_list:
                        self.residues[resnum].shifts[atom] = 0.0
