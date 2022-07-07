#!/usr/bin/python3

"""Validate the "known-good" mmCIF files against the PDBx and MA dictionaries

The files were checked with this script when they were committed, and should
be periodically rechecked in case the PDBx and MA dictionaries are updated.
"""

import ihm.dictionary
import urllib.request

with urllib.request.urlopen(
        'https://raw.githubusercontent.com/ihmwg/ModelCIF/master/'
        'dist/mmcif_ma.dic') as fh:
    pdbx_ma = ihm.dictionary.read(fh)

for to_validate in ('output.cif', 'output_with_align.cif',
                    'test_no_chain.cif', 'test_old_model.cif'):
    with open(to_validate) as fh:
        pdbx_ma.validate(fh)
