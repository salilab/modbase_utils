[![Build Status](https://github.com/salilab/modbase_utils/workflows/build/badge.svg)](https://github.com/salilab/modbase_utils/actions?query=workflow%3Abuild)

This repository contains utilities for interacting with the
[ModBase](https://modbase.compbio.ucsf.edu/) database of
comparative protein structure models:

## `modbase_pdb_to_cif.py`

Convert a PDB file, downloaded from ModBase, to mmCIF format. This should
preserve all information in the PDB file.

This utility requires a local mirror of PDB in compressed mmCIF format.
Use the `-r` command line option to point to the location of this mirror.
It also needs the [python-ihm](https://github.com/ihmwg/python-ihm)
library, to read the mmCIF files for any template structures.
