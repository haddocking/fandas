from fandas.modules.utils import load_bmrbm

from .. import TEST_BMRB_TABLE


def test_load_bmrb():
    """Test the loading of the BMRB table."""
    resnum_col = 1
    atom_col = 3
    shift_col = 4
    observed_dic = load_bmrbm(TEST_BMRB_TABLE, resnum_col, atom_col, shift_col)
    expected_dic = {
        1: {
            "C": 170.52,
            "CA": 54.4,
            "CB": 33.4599,
            "CE": 17.6847,
            "CG": 31.8754,
            "H": 8.2666,
            "HA": 4.22,
            "HB2": 2.07,
            "HB3": 2.139,
            "HE": 1.8684,
            "HG2": 2.4714,
            "HG3": 2.2222,
        },
        2: {"N": 116.08},
        3: {
            "C": 175.22,
            "CA": 54.6,
            "CB": 41.2,
            "CD1": 132.3,
            "CD2": 131.985,
            "CE1": 131.2,
            "CE2": 131.0049,
            "CG": 139.9524,
            "CZ": 129.4,
            "H": 8.57,
        },
        4: {"HG1": 0.71, "HG2": 0.76, "N": 120.7885},
        5: {
            "C": 177.3016,
            "CA": 54.3,
            "CB": 35.5,
            "CD": 29.4536,
            "CE": 42.0471,
            "HG3": 1.4319,
            "HZ": 7.5356,
            "N": 128.38,
        },
    }

    assert expected_dic == observed_dic
