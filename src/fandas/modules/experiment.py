import itertools
import logging

import sys
import copy

from fandas.modules.chemical import ATOM_LIST, ATOM_REF, QUANTUM_RELATIONSHIP

log = logging.getLogger("fandaslog")


class Experiment:
    """Represent the experiments."""

    def __init__(self, chemical_shift, input_data_dic):
        self.chemical_shift = chemical_shift
        self.experiment_dic = input_data_dic["experiments"]
        filtered_shift = None
        distance_file = input_data_dic["distance"]["distance_fname"]
        cutoff = input_data_dic["distance"]["distance_cutoff"]
        if input_data_dic["distance"]["distance_fname"]:
            filtered_shift = self.filter_by_distance(
                distance_fname=distance_file,
                cutoff=cutoff,
            )

        self.filtered_chemical_shift = filtered_shift

    def filter_by_distance(self, distance_fname, cutoff):
        """."""

        distance_dic = self.read_distance_file(distance_fname, cutoff)

        filtered = copy.deepcopy(self.chemical_shift)

        # only keep resnums and atoms that are in contact
        for resnum in self.chemical_shift.residues:
            if resnum in distance_dic.keys():
                for atom in self.chemical_shift.residues[resnum].shifts:
                    if atom not in distance_dic[resnum]:
                        del filtered.residues[resnum].shifts[atom]
        return filtered

    def run_experiments(self):
        """Run all the experiments."""
        for i, exp in enumerate(self.experiment_dic, start=1):

            atoms = self.experiment_dic[exp]["atoms"]
            if "inter" in self.experiment_dic[exp]:
                inter = self.experiment_dic[exp]["inter"]
            else:
                inter = False

            exp_id = "_".join(atoms)

            log.info(f"Running experiment {i} {exp_id}")

            self.execute(atoms, exp_id=exp_id, inter=inter)

            if self.filtered_chemical_shift:
                exp_id += "_dist"
                log.info(f"Running experiment {i} {exp_id}")
                self.execute(atoms, exp_id=exp_id, inter=inter, filtered=True)

    def execute(self, atoms, exp_id, inter=False, filtered=False):
        """Execute the experiment based on an atom list."""
        results = []
        if filtered:
            shifts = self.filtered_chemical_shift
        else:
            shifts = self.chemical_shift

        dimension = len(atoms)
        if "DQ" in atoms:
            if dimension == 2:
                atom_iter = QUANTUM_RELATIONSHIP
        else:
            atom_list = self.translate_atoms(atoms)
            atom_iter = list(itertools.product(*atom_list))

        shift_iterable = self.make_iter(shifts, atom_iter, dimension, inter=inter)

        results = []
        for comb in shift_iterable:
            result_line = self._make_line(comb)
            if result_line:
                print(result_line)
                results.append(result_line)

        if not results:
            log.warning("Nothing to write for this experiment...")
        else:
            self._write_output(results, exp_id)

    def _make_line(self, data_t):
        """."""
        notes = ""
        values = []
        for i, atom_name in zip(data_t[0], data_t[1]):
            notes += f"{i.resname}{i.resnum}{atom_name}-"
            value = self.get_value(i.shifts, atom_name)
            values.append(value)

        if all(values):
            value_str = "\t".join(map("{:.2f}".format, values))
            line = f"{notes[:-1]}\t{value_str}"
            return line

    def make_iter(self, shifts, atom_list, dimension, inter=False):
        """."""
        resnum_list = list(shifts.residues.items())
        first_resum = resnum_list[0][0]
        last_resnum = resnum_list[-1][0]

        # dimension = len(atom_list[0])

        combinations = []
        if dimension == 2:
            for r in range(first_resum, last_resnum + 1):
                r1 = shifts.residues[r]
                if inter:
                    for i in [0, 1, -1]:
                        try:
                            r2 = shifts.residues[r - i]
                            combinations.append((r1, r2))
                        except KeyError:
                            pass
                else:
                    combinations.append((r1, r1))

        elif dimension == 3:
            for r in range(first_resum, last_resnum - 1):
                r1 = shifts.residues[r]
                if inter:
                    for i in [0, 1, -1]:
                        try:
                            r2 = shifts.residues[r - i]
                            combinations.append((r1, r1, r2))
                        except KeyError:
                            pass
                else:
                    combinations.append((r1, r1, r1))

        return list(itertools.product(combinations, atom_list))

    def translate_atoms(self, atoms):
        """Translate the NMR terminology to Python-friendly."""
        translated_atom_l = []
        for raw_atom in atoms:
            if raw_atom in ATOM_REF:
                trans_atom = ATOM_REF[raw_atom]
                if isinstance(trans_atom, str):
                    translated_atom_l.append([trans_atom])
                elif isinstance(trans_atom, list):
                    translated_atom_l.append(trans_atom)
            elif raw_atom in ATOM_LIST:
                translated_atom_l.append([raw_atom])
            else:
                log.error(f"{raw_atom} is unknown.")
                sys.exit()

        return translated_atom_l

    def get_value(self, shift_dic, atom):
        """."""
        return sum([shift_dic[v] for v in atom.split("+")])

    @staticmethod
    def _write_output(results, exp_id):
        """."""
        output_fname = f"{exp_id}_exp_.txt"
        log.info(f"Saving output to {output_fname}")
        with open(output_fname, "w") as fh:
            for line in results:
                fh.write(line)

    @staticmethod
    def read_distance_file(distance_fname, cutoff):
        """."""
        dist_dic = {}
        with open(distance_fname, "r") as fh:
            # 1,CA,1,C,1.548
            for line in fh.readlines():
                data = line.split(",")
                resnum_a, atom_a, resnum_b, atom_b, distance = data
                resnum_a = int(resnum_a)
                # resnum_b = int(resnum_b)
                distance = float(distance)
                if distance <= cutoff:
                    if resnum_a not in dist_dic:
                        dist_dic[resnum_a] = []
                    if atom_a not in dist_dic[resnum_a]:
                        dist_dic[resnum_a].append(atom_a)

        return dist_dic
