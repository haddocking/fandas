"""Test the BMRB module."""
import pytest
from pynmrstar.entry import Entry

from fandas.modules.chemical_shift import BMRB
from fandas.modules.residue import Residue

from .. import TEST_BMRB_TABLE


@pytest.fixture
def bmrb_class():
    """Chemical shift class."""

    yield BMRB(table_fname=TEST_BMRB_TABLE, entity_id=1)


def test__assign(bmrb_class):
    bmrb_class._assign()
    assert bmrb_class.residues[1].shifts == {}
    assert bmrb_class.residues[10].shifts == {"H": 7.82, "N": 109.42}


def test_align_to(bmrb_class):
    aln_dic = bmrb_class.align_to("EPSDTIENVKAKI")
    assert aln_dic == {
        1: 18,
        2: 19,
        3: 20,
        4: 21,
        5: 22,
        6: 23,
        7: 24,
        8: 25,
        9: 26,
        10: 27,
        11: 28,
        12: 29,
        13: 30,
    }


def test__read(bmrb_class):
    entry, residue_dict, seq = bmrb_class._read()

    assert isinstance(entry, Entry)
    assert isinstance(residue_dict, dict)
    assert isinstance(seq, str)
    assert len(residue_dict) == 76
    assert isinstance(residue_dict[1], Residue)
