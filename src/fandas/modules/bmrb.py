import logging
import sys

import pynmrstar
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from pynmrstar.exceptions import ParsingError

from fandas.modules.chemical import THREE_TO_ONE
from fandas.modules.residue import Residue

log = logging.getLogger("fandaslog")


class BMRB:
    def __init__(self, table_fname, entity_id):
        self.table_fname = str(table_fname)
        self.target_entity_id = entity_id
        self.entry, self.residues, self.sequence = self._read()
        self._assign()

    def _assign(self):
        """Assign the BMRB shifts to the residues."""
        for _loop in self.entry.get_loops_by_category("_Atom_chem_shift"):
            for element in _loop.get_tag(
                ["Entity_assembly_ID", "Seq_ID", "Comp_ID", "Atom_ID", "Val"]
            ):
                entity_id, resnum, resname, atom, value = element
                entity_id = int(entity_id)
                if entity_id != self.target_entity_id:
                    continue
                resnum = int(resnum)
                value = float(value)
                resname = THREE_TO_ONE[resname]
                #
                if resnum not in self.residues:
                    self.residues[resnum] = Residue(resnum, resname, None, {})
                self.residues[resnum].shifts[atom] = value

    def align_to(self, reference_sequence):
        """Align the BRMB table to the reference sequence."""
        target_seq = self.sequence
        reference_seq = Seq(reference_sequence)
        aligner = Align.PairwiseAligner(open_gap_score=-10)
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        alns = aligner.align(reference_seq, target_seq)
        top_aln = alns[0]
        identity = (
            str(top_aln).count("|") / float(min(len(reference_seq), len(target_seq)))
        ) * 100
        log.info(f"Identity: {identity:.2f}%")
        for line in str(top_aln).split():
            log.info(line)
        if identity <= 70:
            log.warning("")
            log.warning("!!! Identity is low !!!")
            log.warning("Are you sure this is the correct entity?")
            log.warning("!!! Identity is low !!!")
            log.warning("")

        aligned_ref_segment, aligned_model_segment = top_aln.aligned
        align_dict = {}
        for ref_segment, model_segment in zip(
            aligned_ref_segment, aligned_model_segment
        ):

            start_ref_segment, end_ref_segment = ref_segment
            start_entity_segment, end_entity_segment = model_segment

            for _ref_pos, _entity_pos in zip(
                range(start_ref_segment, end_ref_segment),
                range(start_entity_segment, end_entity_segment),
            ):
                align_dict[_ref_pos + 1] = _entity_pos + 1
        return align_dict

    def _read(self):
        """Use pyNMRStar to read the BMRB table."""
        try:
            entry = pynmrstar.Entry.from_file(self.table_fname)
        except ParsingError as e:
            log.error(e)
            sys.exit()

        residue_dict = {}
        seq = ""
        for _loop in entry.get_loops_by_category("_Entity_poly_seq"):
            for i, element in enumerate(
                _loop.get_tag(["Mon_ID", "Entity_ID"]), start=1
            ):
                resname, entity_id = element
                resname_one_letter = THREE_TO_ONE[resname]
                if int(entity_id) == self.target_entity_id:
                    try:
                        seq += resname_one_letter
                    except KeyError:
                        seq += "X"

                    residue_dict[i] = Residue(i, resname_one_letter, None, {})

        return entry, residue_dict, seq
