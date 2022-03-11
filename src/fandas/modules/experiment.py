from fandas.modules.chemical import (
    atom_type,
    h_atm_ind,
    atom_positions,
    H_ATOMS,
    ATOM_LIST,
)


class Experiment:
    """Represent the experiments."""

    def __init__(self, chemical_shift):
        self.chemical_shift = chemical_shift
        self.sequence = chemical_shift.sequence

    def run(self, name):
        """Run a given experiment."""
        if name == "NH":
            self.nh_exp("H", "N", "H")

    def peaks_proton_heavy(self, target_atom):
        shift_1 = []
        shift_2 = []
        # ALMOST RIGHT, THERE IS AN EXTRA LOOP SOMEWHERE
        #  COMPARE WITH NH_EXP
        for resnum in self.chemical_shift.residues:
            resname = self.chemical_shift.residues[resnum].resname
            for i, atom in enumerate(atom_type):
                if atom == target_atom:
                    for j, h_atom in enumerate(H_ATOMS):
                        # shift of h_atom
                        v1 = self.chemical_shift.residues[resnum].shifts[atom]
                        v2 = self.chemical_shift.residues[resnum].shifts[h_atom]

                        if (
                            # atom == h_atom
                            atom_positions[i] == atom_positions[j]
                            and v1 != 0
                            and v2 != 0
                        ):

                            a = (resname, resnum, v1, ATOM_LIST[j])
                            b = (resname, resnum, v2, ATOM_LIST[i])
                            shift_1.append(a)
                            shift_2.append(b)
        return shift_1, shift_2

    def write_2d(self, shift_a, shift_b, extension, sl=True):
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
        if not peak:
            # log.warning(f"Nothing to write on {extension}_exp.txt")
            pass
        else:
            # log.info(f"Writing {extension}_exp.txt")
            with open(f"{extension}_exp.txt", "w") as output:
                for pk in peak:
                    output.write("%s\t%s\t%s\t%s\n" % (pk[0], pk[1], pk[2], pk[3]))

    def nh_exp(self, atm_1, atm_2, direct):
        shift_1, shift_2 = self.peaks_proton_heavy(atm_2)
        if atm_1 == direct:
            extension = "%s%s" % (atm_1, atm_2)
            self.write_2d(shift_1, shift_2, extension.lower())
        else:
            extension = "%s%s" % (atm_2, atm_1)
            self.write_2d(shift_2, shift_1, extension.lower())
        # ==
        # shift of atom 1

        # pass

        # for i, residue in enumerate(self.sequence):
        #     for j, atom in enumerate(atom_type):
        #         if atom == atm_2:
        #             for h_ind in h_atm_ind:
        #                 atom_positions[h_ind] == atom_positions[j]

        #                 pass
