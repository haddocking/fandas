import logging
from functools import partial

import pandas as pd
from fandas.modules.chemical import (
    AA_REF,
    C_ATOMS,
    DEUTERATION,
    GL_2,
    GL_13,
    N_ATOMS,
    SS_REF,
)
from fandas.modules.residue import Residue
from fandas.modules.utils import load_bmrbm

log = logging.getLogger("fandaslog")


class ChemShift:
    def __init__(self, sequence, secondary_structure):
        self.residues = {}
        self.sequence = sequence
        self.secondary_structure = secondary_structure

    def assign(self, data):
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

    def replace_with_bmrb(
        self, table_fname, resnum_col, atom_col, shift_col, seq_offset
    ):
        """Use a BMRB table to replace the default shift values."""
        bmrb_data = load_bmrbm(table_fname, resnum_col, atom_col, shift_col)

        sequence_number_range = range(seq_offset, len(self.sequence) + 1)
        for bmrb_resnum in bmrb_data:
            # for bmrb_atom in bmrb_data[bmrb_resnum]:
            # atom =
            seq_resnum = (bmrb_resnum + seq_offset) - 1
            seq_resname = self.residues[seq_resnum].resname
            if bmrb_resnum not in sequence_number_range:
                log.warning(
                    f"## bmrb.{bmrb_resnum} does not match sequence numbering, skipping"
                )
                continue

            for bmrb_atom in bmrb_data[bmrb_resnum]:
                new_shift = bmrb_data[bmrb_resnum][bmrb_atom]

                # replace
                try:
                    old_shift = self.residues[seq_resnum].shifts[bmrb_atom]
                    log.info(
                        f"# (bmrb.{bmrb_resnum}) {bmrb_resnum}\t{seq_resname}"
                        f"\t{bmrb_atom}\t{old_shift}\t>> {new_shift}"
                    )

                    self.residues[seq_resnum].shifts[bmrb_atom] = new_shift
                except KeyError:
                    log.warning(
                        f"# (bmrb.{bmrb_resnum}) {bmrb_resnum}\t atom: "
                        f"{bmrb_atom} not found in standard table, skipping..."
                    )

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
        """Apply the foward labeling scheme.

        ...

        """
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
        """Apply the reverse labeling scheme.

        ...

        """
        if rev_12c_14n:
            # remove N and C atoms
            self.zero_shift(C_ATOMS + N_ATOMS, rev_12c_14n)

        if rev_12c:
            self.zero_shift(C_ATOMS, rev_12c)

        if rev_14n:
            self.zero_shift(N_ATOMS, rev_12c)

    def apply_fully_labelling(self):
        """Apply the full labeling."""

        # Do nothing :)

        pass

    def consider_deuteration(self):
        """Consider the deuteration."""
        self.zero_shift(DEUTERATION)

    def zero_shift(self, atoms_to_be_deleted, residue_list=None):
        """Zero the shift of a residue list based on an atom iterable."""
        for resnum in self.residues:
            resname = self.residues[resnum].resname

            if residue_list:
                if resname in residue_list:
                    res_check = True
                else:
                    res_check = False
            else:
                res_check = True

            if isinstance(atoms_to_be_deleted, dict):
                atom_del_list = atoms_to_be_deleted[resname]

            elif isinstance(atoms_to_be_deleted, list):
                atom_del_list = atoms_to_be_deleted

            if res_check:
                for atom in self.residues[resnum].shifts:
                    if atom in atom_del_list:
                        self.residues[resnum].shifts[atom] = 0.0
