#!/usr/bin/python3

"""Validate the "known-good" mmCIF files against the PDBx and MA dictionaries

The files were checked with this script when they were committed, and should
be periodically rechecked in case the PDBx and MA dictionaries are updated.
"""

import ihm.dictionary
import urllib.request

with urllib.request.urlopen(
        'http://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic') as fh:
    pdbx = ihm.dictionary.read(fh)
with urllib.request.urlopen(
        'https://raw.githubusercontent.com/ihmwg/MA-dictionary/'
        'master/mmcif_ma.dic') as fh:
    ma = ihm.dictionary.read(fh)

pdbx_ma = pdbx + ma

for to_validate in ('output.cif', 'output_with_align.cif'):
    with open(to_validate) as fh:
        pdbx_ma.validate(fh)
