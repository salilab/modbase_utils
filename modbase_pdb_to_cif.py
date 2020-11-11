#!/usr/bin/python3

import re


class CifLoop:
    def __init__(self, fh, category, keys):
        self.fh, self.category, self.keys = fh, category, keys
        self._empty_loop = True

    def write(self, line):
        f = self.fh
        if self._empty_loop:
            f.write("#\nloop_\n")
            for k in self.keys:
                f.write("_%s.%s\n" % (self.category, k))
            self._empty_loop = False
        print(line, file=f)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if not self._empty_loop:
            self.fh.write("#\n")


def write_cif_header(fh, model_id, title):
    print("data_model_%s" % model_id, file=fh)
    print("_entry.id modbase_model_%s" % model_id, file=fh)
    print("_struct.entry_id modbase_model_%s" % model_id, file=fh)
    if title:
        print("_struct.title '%s'" % title, file=fh)


def write_cif_exptl(fh, expdta):
    if expdta.startswith('THEORETICAL MODEL, '):
        details = expdta[19:]
        print("#\n_exptl.entry_id modbase_model", file=fh)
        print("_exptl.method 'THEORETICAL MODEL'", file=fh)
        print("_exptl.details '%s'" % details, file=fh)


def write_cif_audit_author(fh):
    with CifLoop(fh, "audit_author", ["name", "pdbx_ordinal"]) as l:
        l.write("'Pieper U' 1'")
        l.write("'Webb B' 2")
        l.write("'Narayanan E' 3")
        l.write("'Sali A' 4")


def write_cif_citation(fh):
    with CifLoop(fh, "citation",
            ["id", "title", "journal_abbrev", "journal_volume", "page_first",
             "page_last", "year", "pdbx_database_id_PubMed",
             "pdbx_database_id_DOI"]) as l:
        l.write("1 'Comparative Protein Modelling by Satisfaction of "
                "Spatial Restraints'\n'J Mol Biol' 234 779 815 1993 8254673 "
                "10.1006/jmbi.1993.1626")
    with CifLoop(fh, "citation_author",
            ["citation_id", "name", "ordinal"]) as l:
        l.write("1 'Sali A' 1")
        l.write("1 'Blundell T L' 2")


def write_cif_software(fh, modpipe_version, modeller_version):
    with CifLoop(fh, "software",
            ["pdbx_ordinal", "name", "classification",
             "version", "type", "location", "citation_id"]) as l:
        l.write("1 ModPipe 'comparative modeling' %s "
                "program https://salilab.org/modpipe/ ." % modpipe_version)
        l.write("2 MODELLER 'comparative modeling' %s program "
                "https://salilab.org/modeller/ 1" % modeller_version)

    # Put each piece of software in its own group
    with CifLoop(fh, "ma_software_group",
            ["ordinal_id", "group_id", "software_id"]) as l:
        for i in range(1,3):
            l.write("%d %d %d" % (i,i,i))


def write_cif_chem_comp(fh):
    # Just assume all 20 standard amino acids are in the model
    with CifLoop(fh, "chem_comp",
            ["id", "type", "name", "formula", "formula_weight"]) as l:
        l.write(
"""ALA 'L-peptide linking' ALANINE 'C3 H7 N O2' 89.094
ARG 'L-peptide linking' ARGININE 'C6 H15 N4 O2 1' 175.212
ASN 'L-peptide linking' ASPARAGINE 'C4 H8 N2 O3' 132.119
ASP 'L-peptide linking' 'ASPARTIC ACID' 'C4 H7 N O4' 133.103
CYS 'L-peptide linking' CYSTEINE 'C3 H7 N O2 S' 121.154
GLN 'L-peptide linking' GLUTAMINE 'C5 H10 N2 O3' 146.146
GLU 'L-peptide linking' 'GLUTAMIC ACID' 'C5 H9 N O4' 147.130
GLY 'peptide linking' GLYCINE 'C2 H5 N O2' 75.067
HIS 'L-peptide linking' HISTIDINE 'C6 H10 N3 O2 1' 156.165
ILE 'L-peptide linking' ISOLEUCINE 'C6 H13 N O2' 131.175
LEU 'L-peptide linking' LEUCINE 'C6 H13 N O2' 131.175
LYS 'L-peptide linking' LYSINE 'C6 H15 N2 O2 1' 147.198
MET 'L-peptide linking' METHIONINE 'C5 H11 N O2 S' 149.208
PHE 'L-peptide linking' PHENYLALANINE 'C9 H11 N O2' 165.192
PRO 'L-peptide linking' PROLINE 'C5 H9 N O2' 115.132
SER 'L-peptide linking' SERINE 'C3 H7 N O3' 105.093
THR 'L-peptide linking' THREONINE 'C4 H9 N O3' 119.120
TRP 'L-peptide linking' TRYPTOPHAN 'C11 H12 N2 O2' 204.229
TYR 'L-peptide linking' TYROSINE 'C9 H11 N O3' 181.191
VAL 'L-peptide linking' VALINE 'C5 H11 N O2' 117.148""")


three_to_one = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LYS' :'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}


def write_cif_entity_details(fh, sequence3):
    # We don't have the alignment, so just one entity (the model, entity_id=1)
    with CifLoop(fh, "entity",
            ["id", "type", "src_method", "pdbx_description"]) as l:
        l.write("1 polymer man target")

    with CifLoop(fh, "entity_poly",
            ["entity_id", "type", "nstd_linkage",
             "pdbx_seq_one_letter_code", "pdbx_seq_one_letter_code_can"]) as l:
        sequence1 = "".join(three_to_one[x] for x in sequence3)
        l.write("1 polypeptide(L) no %s %s" % (sequence1, sequence1))

    with CifLoop(fh, "entity_poly_seq",
            ["entity_id", "num", "mon_id", "hetero"]) as l:
        for i, s in enumerate(sequence3):
            l.write("1 %d %s ." % (i+1, s))


def write_cif_target_details(fh, chain_id):
    with CifLoop(fh, "ma_target_entity",
            ["entity_id", "data_id", "origin"]) as l:
        # target sequence entity=1, data_id=1
        l.write("1 1 .")

    with CifLoop(fh, "ma_target_entity_instance",
            ["asym_id", "entity_id", "details"]) as l:
        l.write("%s 1 ." % chain_id)


def write_cif_assembly(fh, chain_id, sequence3):
    with CifLoop(fh, 'ma_struct_assembly',
            ['ordinal_id', 'assembly_id', 'entity_id', 'asym_id',
             'seq_id_begin', 'seq_id_end']) as l:
        # Simple assembly of a single chain
        l.write("1 1 1 %s 1 %d" % (chain_id, len(sequence3)))


def write_cif_data(fh):
    with CifLoop(fh, "ma_data", ["id", "name", "content_type"]) as l:
        # 1 = target seq, 2 = coordinates
        l.write("1 'Target Sequence' target")
        l.write("2 'Target Structure' 'model coordinates'")

    # Put each data item in its own group
    with CifLoop(fh, "ma_data_group",
            ["ordinal_id", "group_id", "data_id"]) as l:
        for i in range(1, 3):
            l.write("%d %d %d" % (i,i,i))


def write_cif_scores(fh, tsvmod_method, tsvmod_rmsd, tsvmod_no35, ga341, zdope,
        mpqs):
    if not mpqs:
        return
    with CifLoop(fh, 'ma_qa_metric',
            ['id', 'name', 'description', 'type', 'mode',
             'other_details', 'software_group_id']) as l:
        # ModPipe is software_id=1
        l.write("1 MPQS 'ModPipe Quality Score' other global\n"
                "'composite score, values >1.1 are considered reliable' 1")

        # MODELLER is software_id=2
        l.write("2 zDOPE 'Normalized DOPE' zscore global . 2")

        if tsvmod_rmsd:
            l.write("3 'TSVMod RMSD' 'TSVMod predicted RMSD (%s)' "
                    "distance global . ." % tsvmod_method)
            l.write("4 'TSVMod NO35' 'TSVMod predicted native "
                    "overlap (%s)' other global . ." % tsvmod_method)

    with CifLoop(fh, 'ma_qa_metric_global',
            ['ordinal_id', 'model_id', 'metric_id', 'value']) as l:
        l.write("1 1 1 %s" % mpqs)
        l.write("2 1 2 %s" % zdope)
        if tsvmod_rmsd:
            l.write("3 1 3 %s" % tsvmod_rmsd)
            l.write("4 1 4 %s" % tsvmod_no35)


def write_cif_model_list(fh):
    with CifLoop(fh, 'ma_model_list',
            ['ordinal_id', 'model_id', 'model_group_id', 'model_name',
             'model_group_name', 'assembly_id', 'data_id', 'model_type']) as l:
        # Simple output of one model (data_id=2, coordinates)
        l.write("1 1 1 'Selected model' . 1 2 'Homology model'")


def write_cif_asym(fh, chain_id):
    with CifLoop(fh, 'struct_asym', ['id', 'entity_id', 'details']) as l:
        l.write("%s 1 ?" % chain_id)


def write_cif_atom_site(fh, chain_id, atoms, resnum_begin, resnum_end):
    elements = set()
    auth_seqid = resnum_begin
    seqid = 1
    ordinal = 1
    pdb_resnum = None
    with CifLoop(fh, 'atom_site',
            ['group_PDB', 'type_symbol', 'label_atom_id',
             'label_comp_id', 'label_asym_id', 'label_seq_id',
             'auth_seq_id', 'pdbx_PDB_ins_code', 'auth_asym_id',
             'label_alt_id', 'Cartn_x', 'Cartn_y', 'Cartn_z',
             'occupancy', 'B_iso_or_equiv', 'label_entity_id', 'id']) as l:
        for a in atoms:
            # Detect new residue if PDB resnum changed
            pdb_this_resnum = a[22:26]
            if pdb_resnum is not None and pdb_this_resnum != pdb_resnum:
                auth_seqid += 1
                seqid += 1
            pdb_resnum = pdb_this_resnum
            inscode = a[26:27].strip() or '?'
            group_pdb = a[:6]
            element = a[76:78].strip() or '?'
            elements.add(element)
            atmnam = a[12:16]
            resnam = a[17:20]
            x = a[30:38]
            y = a[38:46]
            z = a[46:54]
            occ = a[54:60]
            tfac = a[60:66]
            l.write("%s %s %s %s %s %d %d %s %s . %s %s %s %s %s 1 %d"
                    % (group_pdb, element, atmnam, resnam, chain_id, seqid,
                       auth_seqid, inscode, chain_id, x, y, z, occ, tfac,
                       ordinal))
            ordinal += 1
    assert auth_seqid == resnum_end

    with CifLoop(fh, 'atom_type', ['symbol']) as l:
        for element in sorted(elements):
            l.write(element)


class Structure:
    def _read_pdb(self, fh):
        self.remarks = {}
        self.expdta = None
        self.title = None
        self.modpipe_version = None
        self.atoms = []

        for line in fh:
            if line.startswith('REMARK') and line.count(':') == 1:
                key, val = [x.strip() for x in line[11:].split(':')]
                self.remarks[key] = val
            elif line.startswith('REMARK   6 GENERATED BY MODPIPE VERSION '):
                self.modpipe_version = line[40:].strip()
            elif line.startswith('TITLE     '):
                self.title = line[10:].strip()
            elif line.startswith('EXPDTA    '):
                self.expdta = line[10:].strip()
            elif line.startswith('ATOM') or line.startswith('HETATM'):
                self.atoms.append(line)
        # All ModBase models are single chain
        if self.atoms:
            self.chain_id = self.atoms[0][21]

    def get_modeller_version(self):
        if self.expdta:
            m = re.search('MODELLER\s+(\S+)', self.expdta)
            if m:
                return m.group(1)

    def get_sequence3(self):
        """Get sequence as a sequence of 3-letter residue names"""
        resnum = None
        for a in self.atoms:
            this_resnum = a[22:26]  # residue sequence number
            if this_resnum != resnum:
                yield a[17:20].strip()  # residue name
                resnum = this_resnum

    def write_mmcif(self, fh):
        # mmCIF models must always have a chain ID; older ModBase PDB models had
        # a blank ID
        chain_id = self.chain_id or 'A'
        sequence3 = list(self.get_sequence3())
        modeller_version = self.get_modeller_version() or '?'
        write_cif_header(fh, self.remarks['MODPIPE MODEL ID'], self.title)
        write_cif_exptl(fh, self.expdta)
        write_cif_audit_author(fh)
        write_cif_citation(fh)
        write_cif_software(fh, self.modpipe_version, modeller_version)
        write_cif_chem_comp(fh)
        write_cif_entity_details(fh, sequence3)
        write_cif_target_details(fh, chain_id)
        write_cif_assembly(fh, chain_id, sequence3)
        write_cif_data(fh)
        write_cif_scores(fh, self.remarks.get('TSVMOD METHOD'),
                self.remarks.get('TSVMOD RMSD'),
                self.remarks.get('TSVMOD NO35'),
                self.remarks.get('GA341 SCORE'),
                self.remarks.get('zDOPE SCORE'),
                self.remarks.get('MPQS'))
        write_cif_model_list(fh)
        write_cif_asym(fh, chain_id)
        write_cif_atom_site(fh, chain_id, self.atoms,
                int(self.remarks['TARGET BEGIN']),
                int(self.remarks['TARGET END']))


def read_pdb(fh):
    s = Structure()
    s._read_pdb(fh)
    return s


if __name__ == '__main__':
    import argparse
    a = argparse.ArgumentParser(
            description="Utility to convert ModBase PDB files to mmCIF")
    a.add_argument("pdb", help="Input PDB file")
    a.add_argument("mmcif", help="Output mmCIF file")
    args = a.parse_args()

    with open(args.pdb) as fh:
        s = read_pdb(fh)
    with open(args.mmcif, 'w') as fh:
        s.write_mmcif(fh)
