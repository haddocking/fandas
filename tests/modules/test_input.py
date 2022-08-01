import pytest

from fandas.modules.input import InputFile

from .. import TEST_INPUT_FILE


@pytest.fixture
def inputfile_class():
    """InputFile class."""
    yield InputFile(input_f=TEST_INPUT_FILE)


def test__load(inputfile_class):
    """Test the load method."""
    observed_input = inputfile_class._load(TEST_INPUT_FILE)
    expected_input = {
        "general": {
            "sequence": (
                "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHL"
                "VLRLRGG"
            ),
            "secondary_structure": "",
            "fractionally_deuterated": True,
        },
        "BMRB": {
            "bmrb_table_fname": (""),
            "target_entity": (""),
        },
        "distance": {
            "distance_fname": "",
            "distance_cutoff": 3.9,
        },
        "labeling": {
            "scheme": "fully",
            "forward": {"fw_13c_15n": [], "fw_13c": [], "fw_15n": []},
            "reverse": {"rev_12c_14n": ["A", "G"], "rev_12c": ["L"], "rev_14n": ["I"]},
        },
        "preset_experiments": {
            "selected": [
                "N-H",
                "N-(Calpha)-Cx (residues i, i+1 & i-1)",
                "N-(Co)-Cx",
                "C-C DQ-SQ Correlation",
            ]
        },
        "custom_experiments": {},
    }

    assert observed_input == expected_input


@pytest.mark.skip(reason="Not implemented yet.")
def test__validate(inputfile_class):
    """Test the validate method."""
    inputfile_class._validate()
