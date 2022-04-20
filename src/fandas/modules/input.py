from os import stat
import toml


class InputFile:
    """Represent the toml input file."""

    def __init__(self, input_f):
        self.data = self._load(input_f)

    @staticmethod
    def _validate(data):
        """Validate if the fields are correct."""
        # FIXME: Implement!
        pass

    @staticmethod
    def _load(input_f):
        """Load the input file"""
        return toml.load(input_f)
