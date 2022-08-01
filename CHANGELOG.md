# CHANGELOG

## v2.3.0 (01/08/2022) - Add support to BMRB files in NMR-STAR v3 format

- Only BMRB files in NMR-STAR v3 format are supported.
- If BMRB file is being used, added a feature to automatically align and match the numbering between the input sequence with the sequence observed in the BMRB file.

* * *
v2.2.2 (28/04/2022) - Added missing catalog to the build

* * *
v2.2.1 (26/04/2022) - Bugfix when matching numbering from sequence and chemical shift table

- The matching between the chemical shift data and the sequence numbering was not being done correctly, should not be fixed.

* * *
v2.2.0 (20/04/2022) - Complete code redesign

- Most of the code has been completely redesigned, to make the code more clear, testable and extensible
- Changed the input to a TOML file
- Created a way to represent NMR experiments at `catalog.toml` and changed the code to use this representation
- Previous experiments are now "presets" in the catalog
- Added unittests, coverage report and linting check with github actions

*Developer note: The changes to the code would be worth of a **major** version (v3.0.0) since the API was completely redesigned but in this project it makes more sense for the major version to be incremented when relevant new **scientific** features are added.*

* * *
v2.1.0 (14/02/2022) - Implemented the `offset` feature

- The matching of the numbering between the BMRB table and the sequence was previously done by the server.
- Reformat standard table to a more human-readable format - added `pandas` dependency

* * *
v2.0.3 (12/02/2022) - Fixed bug with indexes in `replace_bmrb`

* * *
v2.0.2 (08/02/2022) - Fixed bug caused by refactoring that caused rev_label function to not recieve the correct arguments

* * *
v2.0.1 (01/02/2022) - Fixed bugs with DQ and distance related experiments
