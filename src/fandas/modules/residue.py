class Residue:
    def __init__(self, resnum, resname, shifts):
        self.resnum = resnum
        self.resname = resname
        self.shifts = shifts

    def __repr__(self):
        return f"Residue({self.resnum} {self.resname})"
