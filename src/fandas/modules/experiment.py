import copy
import itertools
import logging
import os
import sys
from typing import List, Tuple

from fandas.modules.chemical import (
    ATOM_LIST,
    ATOM_REF,
    EXPERIMENT_CATALOG,
    SPECIAL_CASES,
)
from fandas.modules.chemical_shift import ChemShift

log = logging.getLogger("fandaslog")


class Experiment:
    """Represent the experiments."""

    def __init__(self, chemical_shift, input_data_dic):
        self.chemical_shift = chemical_shift
        self.experiment_list = input_data_dic["preset_experiments"]["selected"]
        filtered_shift = None
        distance_file = input_data_dic["distance"]["distance_fname"]
        cutoff = input_data_dic["distance"]["distance_cutoff"]
        self.seq_start = None
        if "BMRB" in input_data_dic.keys():
            if "sequence_start" in input_data_dic["BMRB"].keys():
                self.seq_start = input_data_dic["BMRB"]["sequence_start"]

        if input_data_dic["distance"]["distance_fname"]:
            filtered_shift = self.filter_by_distance(
                distance_fname=distance_file,
                cutoff=cutoff,
            )

        self.filtered_chemical_shift = filtered_shift

    def filter_by_distance(self, distance_fname, cutoff):
        """Filter the chemical shifts by distance."""

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
        for i, nmr_exp_notation in enumerate(self.experiment_list, start=1):
            atoms, direction = self.retrieve_exp_info(nmr_exp_notation)

            log.info(f"Running experiment {i} {nmr_exp_notation}")
            if nmr_exp_notation in SPECIAL_CASES:
                skip_iteration = True
            else:
                skip_iteration = False

            self.execute(
                atoms, direction, exp_id=nmr_exp_notation, skip_iteration=skip_iteration
            )

            if self.filtered_chemical_shift:
                exp_id = f"{nmr_exp_notation}_dist"
                log.info(f"Running experiment {i} {exp_id}")
                self.execute(
                    atoms,
                    direction,
                    exp_id=exp_id,
                    filtered=True,
                    skip_iteration=skip_iteration,
                )

    def execute(self, atoms, direction, exp_id, skip_iteration=False, filtered=False):
        """Execute the experiment based on an atom list."""
        results = []
        if filtered:
            shifts = self.filtered_chemical_shift
        else:
            shifts = self.chemical_shift

        if skip_iteration:
            # Use the pre-defined atom relationships
            atom_iterable = atoms
        else:
            # Generate the atom relationships
            atom_list = self.translate_atoms(atoms)
            atom_iterable = list(itertools.product(*atom_list))

        assert shifts
        shift_iterable = self._make_iter(shifts, atom_iterable, direction)

        results = []
        for comb in shift_iterable:
            result_line = self._make_line(comb)
            if result_line:
                results.append(result_line)

        if not results:
            log.warning("Nothing to write for this experiment...")
        else:
            self._write_output(results, exp_id)

    def _make_line(self, data_t):
        """Generate the line to be written to the output file."""
        notes = ""
        values = []
        # trunk-ignore(ruff/B905)
        for i, atom_name in zip(data_t[0], data_t[1]):
            if self.seq_start is not None:
                notes += f"{i.resname}{self.seq_start + i.resnum - 1}{atom_name}-"
            else:
                notes += f"{i.resname}{i.resnum}{atom_name}-"
            value = self.get_value(i.shifts, atom_name)
            if value != 0.0:
                values.append(value)

        if all(values) and len(values) == len(data_t[0]):
            value_str = "\t".join(map("{:.2f}".format, values))
            line = f"{notes[:-1]}\t{value_str}" + os.linesep

            return line

    def _make_iter(
        self,
        shifts: ChemShift,
        atom_list: List[Tuple[str, str]],
        direction: List[List[int]],
    ) -> List[Tuple[Tuple[ChemShift, ChemShift], str]]:
        """Make an iterable to be used by the experiment executor."""
        resnum_list = list(shifts.residues.items())
        first_resnum = resnum_list[0][0]
        last_resnum = resnum_list[-1][0]

        combinations = []
        for resnum in range(first_resnum, last_resnum + 1):
            for e in itertools.product(*direction):
                # Fetch residues based on adjusted indices
                residues = [
                    shifts.residues[resnum + offset]
                    for offset in e
                    if (resnum + offset) in shifts.residues
                ]

                # Ensure we have the correct number of residues before appending
                if len(residues) == len(e):
                    combinations.append(tuple(residues))
                    # print(e, residues)

        return list(itertools.product(combinations, atom_list))

    def translate_atoms(self, atoms: list):
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
        """Get the shift values."""
        try:
            value = sum([shift_dic[v] for v in atom.split("+")])
        except KeyError:
            value = 0.0

        return value

    @staticmethod
    def retrieve_exp_info(nmr_notation, catalog=EXPERIMENT_CATALOG):
        """Retrieve the atoms and direction from the experiment catalog."""
        atoms = None
        direction = None
        for experiment in catalog:
            if catalog[experiment]["nmr_notation"] == nmr_notation:
                direction = catalog[experiment]["direction"]
                if nmr_notation in SPECIAL_CASES:
                    atoms = catalog[experiment]["atom_relation"]
                else:
                    atoms = catalog[experiment]["atoms"]

        return atoms, direction

    @staticmethod
    def _write_output(results, exp_id):
        """Write the output file."""
        output_fname = f"{exp_id}_exp_.txt"
        log.info(f"Saving output to {output_fname}")
        with open(output_fname, "w") as fh:
            for line in results:
                fh.write(line)

    @staticmethod
    def read_distance_file(distance_fname, cutoff):
        """Read the distance file."""
        dist_dic = {}
        with open(distance_fname, "r") as fh:
            # 1,CA,1,C,1.548
            for line in fh.readlines():
                data = line.split(",")
                resnum_a, atom_a, resnum_b, atom_b, distance = data
                resnum_a = int(resnum_a)
                distance = float(distance)
                if distance <= cutoff:
                    if resnum_a not in dist_dic:
                        dist_dic[resnum_a] = []
                    if atom_a not in dist_dic[resnum_a]:
                        dist_dic[resnum_a].append(atom_a)

        return dist_dic
