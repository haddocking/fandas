"""Test module."""
from pathlib import Path

TEST_DATA_PATH = Path(Path(__file__).resolve().parents[0], "data")
TEST_BMRB_TABLE = Path(TEST_DATA_PATH, "test_bmrm_table.csv")
TEST_INPUT_FILE = Path(TEST_DATA_PATH, "test_input.toml")
TEST_DISTANCE_FILE = Path(TEST_DATA_PATH, "test_distance.txt")
