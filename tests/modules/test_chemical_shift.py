"""Test the chemical shift module."""

import copy

import pytest
from fandas.modules.chemical_shift import ChemShift

from .. import TEST_BMRB_TABLE


@pytest.fixture
def chemshift_class():
    """Chemical shift class."""
    sequence = "MQIFV"
    secondary_structure = "naabc"

    yield ChemShift(sequence, secondary_structure)


def test_assign(chemshift_class):
    """Test the assign method."""
    chemshift_class.assign()

    assert len(chemshift_class.residues) == 5

    expected_ss = [
        "COIL",
        "ALPHA",
        "ALPHA",
        "BETA",
        "COIL",
    ]

    for i, resnum in enumerate(chemshift_class.residues):
        assert expected_ss[i] == chemshift_class.residues[resnum].secondary_structure

    assert chemshift_class.residues[1].resname == "M"
    assert chemshift_class.residues[1].shifts["H"] == 8.37
    assert chemshift_class.residues[1].shifts["CE"] == 17.25
    assert chemshift_class.residues[1].shifts["C"] == 175.93


def test_replace_with_bmrb(chemshift_class):
    """Test the replace_with_bmrb method."""
    chemshift_class.replace_with_bmrb(
        table_fname=TEST_BMRB_TABLE, resnum_col=1, atom_col=3, shift_col=4, seq_offset=1
    )

    assert chemshift_class.residues[1].shifts["C"] == 170.52
    assert chemshift_class.residues[1].shifts["CA"] == 54.4
    assert chemshift_class.residues[5].shifts["N"] == 128.38


@pytest.mark.skip(reason="No way of currently testing this")
def test_label(chemshift_class):
    pass


@pytest.mark.skip(reason="No way of currently testing this")
def test__get_label_func(chemshift_class):
    pass


def test_apply_gl13_labelling(chemshift_class):
    """Test the apply_gl13_labelling method."""
    chemshift_class.apply_gl13_labelling()
    assert chemshift_class.residues[5].shifts["CA"] == 0.0


def test_apply_gl2_labelling(chemshift_class):
    """Test the apply_gl2_labelling method."""
    chemshift_class.apply_gl2_labelling()
    assert chemshift_class.residues[5].shifts["C"] == 0.0


def test_apply_fw_labelling(chemshift_class):
    """Test the apply_fw_labelling method."""
    chemshift_class.apply_fw_labelling(fw_13c_15n=[1, 2], fw_13c=[3], fw_15n=[4])

    assert chemshift_class.residues[1].shifts["CA"] == 0.0
    assert chemshift_class.residues[2].shifts["CA"] == 0.0

    assert chemshift_class.residues[1].shifts["C"] == 0.0
    assert chemshift_class.residues[1].shifts["N"] == 0.0

    assert chemshift_class.residues[3].shifts["N"] == 120.22
    assert chemshift_class.residues[3].shifts["CA"] == 0.0

    assert chemshift_class.residues[4].shifts["CA"] == 56.33
    assert chemshift_class.residues[4].shifts["N"] == 0.0


def test_apply_rev_labelling(chemshift_class):
    """Test the apply_rev_labelling method."""
    chemshift_class.apply_rev_labelling(rev_12c_14n=[1, 2], rev_12c=[3], rev_14n=[4])

    assert chemshift_class.residues[1].shifts["CA"] == 0.0
    assert chemshift_class.residues[2].shifts["CA"] == 0.0

    assert chemshift_class.residues[1].shifts["C"] == 0.0
    assert chemshift_class.residues[1].shifts["N"] == 0.0

    assert chemshift_class.residues[3].shifts["N"] == 120.22
    assert chemshift_class.residues[3].shifts["CA"] == 0.0

    assert chemshift_class.residues[4].shifts["CA"] == 56.33
    assert chemshift_class.residues[4].shifts["N"] == 0.0


def test_apply_fully_labelling(chemshift_class):
    """Test the apply_fully_labelling method."""
    unmodified = copy.deepcopy(chemshift_class)

    chemshift_class.apply_fully_labelling()

    assert unmodified.residues[1].shifts == chemshift_class.residues[1].shifts
    assert unmodified.residues[2].shifts == chemshift_class.residues[2].shifts
    assert unmodified.residues[3].shifts == chemshift_class.residues[3].shifts
    assert unmodified.residues[4].shifts == chemshift_class.residues[4].shifts
    assert unmodified.residues[5].shifts == chemshift_class.residues[5].shifts


def test_consider_deuteration(chemshift_class):
    """Test the consider_deuteration method."""
    chemshift_class.consider_deuteration()

    assert chemshift_class.residues[1].shifts["HA"] == 0.0
    assert chemshift_class.residues[2].shifts["HA"] == 0.0
    assert chemshift_class.residues[3].shifts["HB"] == 0.0
    assert chemshift_class.residues[3].shifts["HB"] == 0.0


def test_zero_shift(chemshift_class):
    """Test the zero_shift method."""
    chemshift_class.zero_shift(atoms_to_be_set_to_zero="C")
    assert chemshift_class.residues[1].shifts["C"] == 0.0

    chemshift_class.zero_shift(atoms_to_be_set_to_zero=["CA"])
    assert chemshift_class.residues[1].shifts["CA"] == 0.0

    chemshift_class.zero_shift(atoms_to_be_set_to_zero={"I": ["CG2"]})
    assert chemshift_class.residues[3].shifts["CG2"] == 0.0
