import os
import unittest
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
INPUTDIR = os.path.join(TOPDIR, 'test', 'input')

inpdb = os.path.join(INPUTDIR,
        'modbase-model_66b8fbc891f519c1ba8d8ad2e62c6caa.pdb')
align = os.path.join(INPUTDIR, 'align.ali')
script = os.path.join(TOPDIR, 'modbase_pdb_to_cif.py')

# Known-good output mmCIF files. These should be periodically checked
# for validity using input/validate-modbase.py
outcif = os.path.join(INPUTDIR, 'output.cif')
outaligncif = os.path.join(INPUTDIR, 'output_with_align.cif')


class Tests(unittest.TestCase):
    def compare_files(self, f1, f2):
        with open(f1) as fh:
            f1_contents = fh.read()
        with open(f2) as fh:
            f2_contents = fh.read()
        self.assertEqual(f1_contents, f2_contents)

    def test_without_align(self):
        """Test modbase_pdb_to_cif script without alignment"""
        out = 'test_without_align.cif'
        p = subprocess.check_call([script, inpdb, out])
        self.compare_files(out, outcif)
        os.unlink(out)

    def test_with_align(self):
        """Test modbase_pdb_to_cif script with alignment"""
        out = 'test_with_align.cif'
        p = subprocess.check_call([script, '-a', align, inpdb, out])
        self.compare_files(out, outaligncif)
        os.unlink(out)


if __name__ == '__main__':
    unittest.main()
