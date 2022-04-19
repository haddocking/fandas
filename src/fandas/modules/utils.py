import logging

log = logging.getLogger("fandaslog")


def load_bmrbm(table_fname, resname_col, atom_col, shift_col):
    """Load a BMRBM table and return a dictionary."""
    bmrb_dic = {}
    with open(table_fname, "r") as bmrb_file:
        for line in bmrb_file.readlines():

            # Split so that the whitespaces don't matter
            c_shift = line.split()

            if not c_shift:
                # line is empty
                continue

            #  bmrb_columns starts with index 1 so here we need to -1 all positions
            bmrb_resnum = int(c_shift[resname_col - 1])
            bmrb_atom = c_shift[atom_col - 1]
            bmrb_shift = float(c_shift[shift_col - 1])

            if bmrb_resnum not in bmrb_dic:
                bmrb_dic[bmrb_resnum] = {}

            bmrb_dic[bmrb_resnum][bmrb_atom] = bmrb_shift

    return bmrb_dic
