import toml


class InputFile:
    """Represent the toml input file."""

    def __init__(self, input_f):
        self.data = self._load(input_f)

    def _load(self, input_f):
        """Load the input file"""
        data = toml.load(input_f)
        return self._validate(data)

    @staticmethod
    def _validate(data):
        """Validate if the fields are correct."""
        # FIXME: Implement!
        return data
