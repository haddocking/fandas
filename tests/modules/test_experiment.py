"""Test the experiment module."""
import pytest
import tempfile
import os
from pathlib import Path
from fandas.modules.experiment import Experiment
from fandas.modules.chemical_shift import ChemShift
from fandas.modules.input import InputFile

from .. import TEST_INPUT_FILE, TEST_DISTANCE_FILE


@pytest.fixture
def chemical_shifts():
    """Chemical shift class."""
    sequence = "MQIFV"
    secondary_structure = "naabc"

    yield ChemShift(sequence, secondary_structure)


@pytest.fixture
def inputfile_data():
    """InputFile class."""
    yield InputFile(input_f=TEST_INPUT_FILE).data


@pytest.fixture
def experiment_class(chemical_shifts, inputfile_data):
    """Experiment class."""
    yield Experiment(chemical_shift=chemical_shifts, input_data_dic=inputfile_data)


def test_filter_by_distance(experiment_class):
    """Test the filter_by_distance method."""
    filtered = experiment_class.filter_by_distance(
        distance_fname=TEST_DISTANCE_FILE, cutoff=2.0
    )

    assert "C" not in filtered.residues[1].shifts
    assert "CA" in filtered.residues[1].shifts


def test_run_experiments(experiment_class):
    """Test the run_experiments method."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        os.chdir(tmpdirname)

        experiment_class.run_experiments()

        # TODO: run ALL experiments here
        assert Path("C-C DQ-SQ Correlation_exp_.txt").exists()
        assert Path("N-(Calpha)-Cx (residues i, i+1 & i-1)_exp_.txt").exists()
        assert Path("N-H_exp_.txt").exists()
        assert Path("N-(Co)-Cx_exp_.txt").exists()


def test_execute(experiment_class):
    """Test the execute method."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        os.chdir(tmpdirname)
        experiment_class.execute(
            atoms=["H", "H"], direction=[[0], [0]], exp_id="test_2d"
        )
        assert Path("test_2d_exp_.txt").exists()
        assert Path("test_2d_exp_.txt").stat().st_size > 0

        experiment_class.execute(
            atoms=["H", "H", "H"], direction=[[0], [0], [-1, 0]], exp_id="test_3d"
        )
        assert Path("test_3d_exp_.txt").exists()
        assert Path("test_3d_exp_.txt").stat().st_size > 0


def test__make_line(experiment_class, chemical_shifts):
    """Test the _make_line method."""
    res_a = chemical_shifts.residues[1]
    res_b = chemical_shifts.residues[1]
    observed_line = experiment_class._make_line(data_t=((res_a, res_b), ("N", "H")))
    expected_line = "M1N-M1H\t120.19\t8.37" + os.linesep

    assert observed_line == expected_line


def test__make_iter(experiment_class, chemical_shifts):
    """Test the _make_iter method."""
    observed_iter = experiment_class._make_iter(
        shifts=chemical_shifts, atom_list=[("H", "H")], direction=[[0], [0]]
    )
    expected_iter = [
        ((chemical_shifts.residues[1], chemical_shifts.residues[1]), ("H", "H")),
        ((chemical_shifts.residues[2], chemical_shifts.residues[2]), ("H", "H")),
        ((chemical_shifts.residues[3], chemical_shifts.residues[3]), ("H", "H")),
        ((chemical_shifts.residues[4], chemical_shifts.residues[4]), ("H", "H")),
        ((chemical_shifts.residues[5], chemical_shifts.residues[5]), ("H", "H")),
    ]

    assert observed_iter == expected_iter

    observed_iter = experiment_class._make_iter(
        shifts=chemical_shifts, atom_list=[("H", "H"), ("H", "N")], direction=[[0], [0]]
    )
    expected_iter = [
        ((chemical_shifts.residues[1], chemical_shifts.residues[1]), ("H", "H")),
        ((chemical_shifts.residues[1], chemical_shifts.residues[1]), ("H", "N")),
        ((chemical_shifts.residues[2], chemical_shifts.residues[2]), ("H", "H")),
        ((chemical_shifts.residues[2], chemical_shifts.residues[2]), ("H", "N")),
        ((chemical_shifts.residues[3], chemical_shifts.residues[3]), ("H", "H")),
        ((chemical_shifts.residues[3], chemical_shifts.residues[3]), ("H", "N")),
        ((chemical_shifts.residues[4], chemical_shifts.residues[4]), ("H", "H")),
        ((chemical_shifts.residues[4], chemical_shifts.residues[4]), ("H", "N")),
        ((chemical_shifts.residues[5], chemical_shifts.residues[5]), ("H", "H")),
        ((chemical_shifts.residues[5], chemical_shifts.residues[5]), ("H", "N")),
    ]

    observed_iter = experiment_class._make_iter(
        shifts=chemical_shifts,
        atom_list=[("H", "H", "H"), ("H", "N", "H")],
        direction=[[0], [0], [+1, 0]],
    )

    expected_iter = [
        (
            (
                chemical_shifts.residues[1],
                chemical_shifts.residues[1],
                chemical_shifts.residues[1],
            ),
            ("H", "H", "H"),
        ),
        (
            (
                chemical_shifts.residues[1],
                chemical_shifts.residues[1],
                chemical_shifts.residues[1],
            ),
            ("H", "N", "H"),
        ),
        (
            (
                chemical_shifts.residues[2],
                chemical_shifts.residues[2],
                chemical_shifts.residues[1],
            ),
            ("H", "H", "H"),
        ),
        (
            (
                chemical_shifts.residues[2],
                chemical_shifts.residues[2],
                chemical_shifts.residues[1],
            ),
            ("H", "N", "H"),
        ),
        (
            (
                chemical_shifts.residues[2],
                chemical_shifts.residues[2],
                chemical_shifts.residues[2],
            ),
            ("H", "H", "H"),
        ),
        (
            (
                chemical_shifts.residues[2],
                chemical_shifts.residues[2],
                chemical_shifts.residues[2],
            ),
            ("H", "N", "H"),
        ),
        (
            (
                chemical_shifts.residues[3],
                chemical_shifts.residues[3],
                chemical_shifts.residues[2],
            ),
            ("H", "H", "H"),
        ),
        (
            (
                chemical_shifts.residues[3],
                chemical_shifts.residues[3],
                chemical_shifts.residues[2],
            ),
            ("H", "N", "H"),
        ),
        (
            (
                chemical_shifts.residues[3],
                chemical_shifts.residues[3],
                chemical_shifts.residues[3],
            ),
            ("H", "H", "H"),
        ),
        (
            (
                chemical_shifts.residues[3],
                chemical_shifts.residues[3],
                chemical_shifts.residues[3],
            ),
            ("H", "N", "H"),
        ),
        (
            (
                chemical_shifts.residues[4],
                chemical_shifts.residues[4],
                chemical_shifts.residues[3],
            ),
            ("H", "H", "H"),
        ),
        (
            (
                chemical_shifts.residues[4],
                chemical_shifts.residues[4],
                chemical_shifts.residues[3],
            ),
            ("H", "N", "H"),
        ),
        (
            (
                chemical_shifts.residues[4],
                chemical_shifts.residues[4],
                chemical_shifts.residues[4],
            ),
            ("H", "H", "H"),
        ),
        (
            (
                chemical_shifts.residues[4],
                chemical_shifts.residues[4],
                chemical_shifts.residues[4],
            ),
            ("H", "N", "H"),
        ),
        (
            (
                chemical_shifts.residues[5],
                chemical_shifts.residues[5],
                chemical_shifts.residues[4],
            ),
            ("H", "H", "H"),
        ),
        (
            (
                chemical_shifts.residues[5],
                chemical_shifts.residues[5],
                chemical_shifts.residues[4],
            ),
            ("H", "N", "H"),
        ),
        (
            (
                chemical_shifts.residues[5],
                chemical_shifts.residues[5],
                chemical_shifts.residues[5],
            ),
            ("H", "H", "H"),
        ),
        (
            (
                chemical_shifts.residues[5],
                chemical_shifts.residues[5],
                chemical_shifts.residues[5],
            ),
            ("H", "N", "H"),
        ),
    ]

    assert observed_iter == expected_iter


def test_translate_atoms(experiment_class):
    """Test the translation of atoms."""
    observed_translation = experiment_class.translate_atoms(["H", "CO"])
    expected_translation = [["H"], ["C"]]
    assert observed_translation == expected_translation


def test_get_value(experiment_class):
    """Test the retrieval of values."""
    observed_value = experiment_class.get_value(shift_dic={"H": 1.0}, atom="H")
    expected_value = 1.0

    assert observed_value == expected_value

    observed_value = experiment_class.get_value(
        shift_dic={"H": 1.0, "N": 2.0}, atom="H+N"
    )
    expected_value = 3.0

    assert observed_value == expected_value


def test_retrieve_exp_info(experiment_class):
    """Test the retrieval of experimental information."""
    observed_atoms, observed_direction = experiment_class.retrieve_exp_info(
        nmr_notation="N-H"
    )
    expected_atoms = ["N", "H"]
    expected_direction = [[0], [0]]

    assert observed_atoms == expected_atoms
    assert observed_direction == expected_direction

    observed_atoms, observed_direction = experiment_class.retrieve_exp_info(
        nmr_notation="C-C DQ-SQ Correlation"
    )

    expected_direction = [[0], [0]]

    # FIXME: Add a proper check for the quantum relatiom here
    assert observed_atoms is not None
    assert observed_direction == expected_direction


def test__write_output(experiment_class):
    """Test the writing of output."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        os.chdir(tmpdirname)
        experiment_class._write_output(results=["TestString"], exp_id="test")

        assert Path("test_exp_.txt").exists()
        assert Path("test_exp_.txt").read_text() == "TestString"
        assert Path("test_exp_.txt").stat().st_size > 0


def test_read_distance_file(experiment_class):
    """Test the reading of distance files."""
    observed_distance_dic = experiment_class.read_distance_file(
        distance_fname=TEST_DISTANCE_FILE, cutoff=2.0
    )
    expected_distance_dic = {1: ["CA"]}

    assert observed_distance_dic == expected_distance_dic
