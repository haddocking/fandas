import toml


class InputFile:
    """Represent the toml input file."""

    def __init__(self, input_f):
        self.data = self._load(input_f)
        self._validate()

    def _validate(self):
        """Validate if the fields are correct."""
        # FIXME: Implement!

        # validate BMRB fields
        if "BMRB" not in self.data:
            # these can be empty
            self.data["BMRB"] = {}
            self.data["BMRB"]["bmrb_table_fname"] = None
            self.data["BMRB"]["target_entity"] = None
        else:
            if "bmrb_table_fname" not in self.data["BMRB"]:
                raise ValueError("Missing `bmrb_table_fname` parameter in BMRB section")
            if "target_entity" not in self.data["BMRB"]:
                raise ValueError("Missing `target_entity` parameter in BMRB section")

    @staticmethod
    def _load(input_f):
        """Load the input file"""
        return toml.load(input_f)
