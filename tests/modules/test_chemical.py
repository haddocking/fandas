"""Test the chemical module."""
from fandas.modules.chemical import (
    STANDARD_DATA,
    ATOM_LIST,
    QUANTUM_RELATIONSHIP,
    EXPERIMENT_CATALOG,
    ATOM_INDEX_DIC,
    N_ATOMS,
    C_ATOMS,
    GL_13,
    AA_REF,
    SS_REF,
    GL_2,
    DEUTERATION,
)


def test_standard_data():
    """Test if the standard data is defines."""
    assert STANDARD_DATA is not None


def test_atom_list():
    """Test if the atom list is defined."""
    assert ATOM_LIST is not None


def test_quantum_relationship():
    """Test if the quantum relationship is defined."""
    assert QUANTUM_RELATIONSHIP is not None


def test_experiment_catalog():
    """Test if the experiment catalog is defined."""
    assert EXPERIMENT_CATALOG is not None

    # TODO: add a check for each experiment


def test_atom_index_dic():
    """Test if the atom index dictionary is defined."""
    assert ATOM_INDEX_DIC is not None


def test_n_atoms():
    """Test if the N_ATOMS is defined."""
    assert N_ATOMS is not None
    assert len(N_ATOMS) == 9


def test_c_atoms():
    """Test if the C_ATOMS is defined."""
    assert C_ATOMS is not None
    assert len(C_ATOMS) == 17


def test_gl13():
    """Test if the GL13 atom dictionary is defined."""
    assert GL_13 is not None


def test_aa_ref():
    """Test if the AA_REF is defined."""
    assert AA_REF is not None


def test_ss_ref():
    """Test if the SS_REF is defined."""
    assert SS_REF is not None


def test_gl2():
    """Test if the GL_REF is defined."""
    assert GL_2 is not None


def test_deuteration():
    """Test if the DEUTERATION is defined."""
    assert DEUTERATION is not None
