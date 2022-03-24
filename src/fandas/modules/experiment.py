import itertools
import logging
import os
import sys

from fandas.modules.chemical import ATOM_LIST, ATOM_REF

log = logging.getLogger("fandaslog")


class Experiment:
    """Represent the experiments."""

    def __init__(self, chemical_shift, experiment_dic):
        self.chemical_shift = chemical_shift
        self.experiment_dic = experiment_dic

    def run(self):
        """Run all the experiments."""
        for i, exp in enumerate(self.experiment_dic, start=1):

            atoms = self.experiment_dic[exp]
            exp_id = "_".join(atoms)

            log.info(f"Running experiment {i} {exp_id}")

            atoms = self.translate_atoms(atoms)

            self.execute(atoms, exp_id=exp_id)

    def execute(self, atom_list, exp_id):
        """Execute the experiment based on an atom list."""
        results = []
        for resnum in self.chemical_shift.residues:
            residue = self.chemical_shift.residues[resnum]
            try:
                value_list = [residue.shifts[atom] for atom in atom_list]
                if all(value_list):
                    result_line = self._make_line(residue, atom_list, value_list)
                    results.append(result_line)
            except TypeError:
                # this is a one-to-many experiment
                # N-C, N-CB, N-CG, etc
                for comb in itertools.product(*atom_list):
                    value_list = [residue.shifts[atom] for atom in comb]
                    if all(value_list):
                        result_line = self._make_line(residue, comb, value_list)
                        results.append(result_line)

        if not results:
            log.warning("Nothing to write for this experiment...")
        else:
            self._write_output(results, exp_id)

    def translate_atoms(self, atoms):
        """Translate the NMR terminology to Python-friendly."""
        translated_atom_l = []
        for raw_atom in atoms:
            if raw_atom in ATOM_REF:
                trans_atom = ATOM_REF[raw_atom]
                translated_atom_l.append(trans_atom)
            elif raw_atom in ATOM_LIST:
                translated_atom_l.append(raw_atom)
            else:
                log.error(f"{raw_atom} is unknown.")
                sys.exit()

        return translated_atom_l

    @staticmethod
    def _write_output(results, exp_id):
        """."""
        output_fname = f"{exp_id}_exp_.txt"
        log.info(f"Saving output to {output_fname}")
        with open(output_fname, "w") as fh:
            for line in results:
                fh.write(line)

    @staticmethod
    def _make_line(residue, atom_list, value_list):
        """."""
        identifier = "-".join(
            [f"{residue.resname}{residue.resnum}{atom}" for atom in atom_list]
        )
        values_str = "\t".join(map(str, value_list))
        return f"{identifier}\t{values_str}" + os.linesep
