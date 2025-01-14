import os
import sys
import unittest
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
INPUTDIR = os.path.join(TOPDIR, 'test', 'input')

inpdb = os.path.join(
    INPUTDIR, 'modbase-model_66b8fbc891f519c1ba8d8ad2e62c6caa.pdb')
align = os.path.join(INPUTDIR, 'align.ali')
repo = os.path.join(INPUTDIR, 'mmcif_repo')
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
        _ = subprocess.check_call([sys.executable, script, '-r', repo,
                                   inpdb, out])
        self.compare_files(out, outcif)
        os.unlink(out)

    def test_with_align(self):
        """Test modbase_pdb_to_cif script with alignment"""
        out = 'test_with_align.cif'
        _ = subprocess.check_call([sys.executable, script, '-r', repo,
                                   '-a', align, inpdb, out])
        self.compare_files(out, outaligncif)
        os.unlink(out)

    def test_with_align_missing_template(self):
        """Test modbase_pdb_to_cif with alignment but missing template"""
        out = 'test_with_align_missing_template.cif'
        _ = subprocess.check_call([sys.executable, script,
                                   '-a', align, inpdb, out])
        # Without the repo, we won't be able to find the template .cif;
        # information on asyms and residue ranges will be missing
        with open(out) as fh:
            f1_contents = fh.read()
        with open(outaligncif) as fh:
            f2_contents = fh.read()
        f2_contents = f2_contents.replace('polymer 1 1 A A 1 1 T',
                                          'polymer 1 1 A ? ? 1 T')
        f2_contents = f2_contents.replace('1 1 403 837', '1 1 ? ?')
        self.assertEqual(f1_contents, f2_contents)
        os.unlink(out)

    def test_with_align_bad_range(self):
        """Test modbase_pdb_to_cif with alignment and bad sequence range"""
        out = 'test_with_align_bad_range.cif'
        # Modify input PDB to reference a bad template sequence range
        inpdb_bad_range = 'inpdb_bad_range.pdb'
        with open(inpdb) as fh:
            contents = fh.read()
        contents = contents.replace('TEMPLATE BEGIN:              401',
                                    'TEMPLATE BEGIN:              901')
        with open(inpdb_bad_range, 'w') as fh:
            fh.write(contents)
        _ = subprocess.check_call([sys.executable, script, '-r', repo,
                                   '-a', align, inpdb_bad_range, out])
        # Information on asyms and residue ranges will be missing
        with open(out) as fh:
            f1_contents = fh.read()
        with open(outaligncif) as fh:
            f2_contents = fh.read()
        f2_contents = f2_contents.replace('polymer 1 1 A A 1 1 T',
                                          'polymer 1 1 A ? ? 1 T')
        f2_contents = f2_contents.replace('1 1 403 837', '1 1 ? ?')
        self.assertEqual(f1_contents, f2_contents)

        os.unlink(inpdb_bad_range)
        os.unlink(out)

    def test_bulk_headers(self):
        """Test modbase_pdb_to_cif script without bulk-headers PDB"""
        inpdb = os.path.join(INPUTDIR, 'test_bulk_headers.pdb')
        out = 'test_without_align.cif'
        _ = subprocess.check_call([sys.executable, script, '-r', repo,
                                   inpdb, out])
        self.compare_files(out, outcif)
        os.unlink(out)

    def test_no_chain(self):
        """Test modbase_pdb_to_cif script with a PDB with no chain ID"""
        out = 'test_no_chain.cif'
        inpdb = os.path.join(INPUTDIR, 'test_no_chain.pdb')
        outcif = os.path.join(INPUTDIR, 'test_no_chain.cif')
        _ = subprocess.check_call([sys.executable, script, '-r', repo,
                                   inpdb, out])
        self.compare_files(out, outcif)
        os.unlink(out)

    def test_old_model(self):
        """Test modbase_pdb_to_cif script with old-style models"""
        out = 'test_old_model.cif'
        inpdb = os.path.join(INPUTDIR, 'test_old_model.pdb')
        align = os.path.join(INPUTDIR, 'test_old_model.ali')
        outcif = os.path.join(INPUTDIR, 'test_old_model.cif')
        _ = subprocess.check_call([sys.executable, script, '-r', repo,
                                   '-a', align, inpdb, out])
        self.compare_files(out, outcif)
        os.unlink(out)

    def test_multi_line_and_expdta(self):
        """Test with multiline TITLE and new-style EXPDTA record"""
        out = 'test_multi_expdta.cif'
        inpdb = os.path.join(INPUTDIR, 'test_multi_expdta.pdb')
        outcif = os.path.join(INPUTDIR, 'test_multi_expdta.cif')
        _ = subprocess.check_call([sys.executable, script, '-r', repo,
                                   inpdb, out])
        self.compare_files(out, outcif)
        os.unlink(out)


if __name__ == '__main__':
    unittest.main()
