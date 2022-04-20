import pytest
from fandas.modules.residue import Residue


@pytest.fixture
def residue_class():
    """Residue class."""
    yield Residue(
        resnum=1, resname="M", secondary_structure="ALPHA", shifts={"CA": 54.4}
    )


def test___repr___(residue_class):
    """Test the __repr__ method."""
    assert str(residue_class) == "Residue(1.M.ALPHA)"
