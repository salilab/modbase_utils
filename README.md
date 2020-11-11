This repository contains utilities for interacting with the
[ModBase](https://modbase.compbio.ucsf.edu/) database of
comparative protein structure models:

## `modbase_pdb_to_cif.py`

Convert a PDB file, downloaded from ModBase, to mmCIF format. This should
preserve all information in the PDB file.

Note that the resulting mmCIF file differs from an mmCIF downloaded directly
from ModBase in that it contains no information on the modeling alignment
(since this information is not present in ModBase PDB files).
