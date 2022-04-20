class Residue:
    def __init__(self, resnum, resname, secondary_structure, shifts):
        self.resnum = resnum
        self.resname = resname
        self.secondary_structure = secondary_structure
        self.shifts = shifts

    def __repr__(self):
        return f"Residue({self.resnum}.{self.resname}.{self.secondary_structure})"
