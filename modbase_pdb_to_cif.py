#!/usr/bin/python3

import re
import collections

# Single sequence in a Modeller alignment
Sequence = collections.namedtuple("Sequence",
        ["seqtyp", "chain", "method", "gapped", "primary"])

# Mapping between one-letter codes and PDB names
three_to_one = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LYS' :'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}
one_to_three = {val:key for key, val in three_to_one.items()}


class Alignment:
    """Represent a Modeller alignment"""

    def __init__(self, fname):
        with open(fname) as fh:
            self.template, self.target = self._read_seqs(fh)

    def _read_seqs(self, fh):
        template, target = None, None
        for line in fh:
            if line.startswith('>P1;'):
                seq = self._read_seq(fh)
                if seq.seqtyp == 'sequence':
                    target = seq
                elif seq.seqtyp.startswith('structure'):
                    template = seq
        if template is None or target is None:
            raise ValueError("Could not read target and template "
                             "from alignment")
        # All current ModBase models have only template
        return template, target

    def _read_seq(self, fh):
        header = fh.readline().split(':')
        seqlines = []
        while True:
            line = fh.readline()
            if line == '':
                raise ValueError("End of file while reading sequence")
            seqlines.append(line.rstrip('\r\n'))
            if seqlines[-1].endswith('*'):
                break
        gapped = "".join(seqlines)[:-1]
        return Sequence(seqtyp=header[0], chain=header[3], method=header[7],
                gapped=gapped, primary=gapped.replace('-', ''))


class CifLoop:
    """Helper class to write an mmCIF loop construct"""

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


class CifWriter:
    def __init__(self, fh, align):
        class CifID:
            pass

        self.fh, self.align = fh, align
        # Assign consecutive CIF numeric IDs from 1
        self.target = CifID()
        self.coord = CifID()
        self.template = CifID()
        self.alignment = CifID()
        if align:
            self.template.entity_id, self.target.entity_id = 1, 2
        else:
            self.target.entity_id = 1
        self.template.data_id, self.target.data_id = 1, 2
        self.alignment.data_id, self.coord.data_id = 3, 4

    def print(self, s):
        print(s, file=self.fh)

    def loop(self, category, keys):
        return CifLoop(self.fh, category, keys)

    def write_header(self, model_id, title):
        self.print("data_model_%s" % model_id)
        self.print("_entry.id modbase_model_%s" % model_id)
        self.print("_struct.entry_id modbase_model_%s" % model_id)
        if title:
            self.print("_struct.title '%s'" % title)

    def write_exptl(self, model_id, expdta):
        if expdta.startswith('THEORETICAL MODEL, '):
            details = expdta[19:]
            self.print("#\n_exptl.entry_id modbase_model_%s" % model_id)
            self.print("_exptl.method 'THEORETICAL MODEL'")
            self.print("_exptl.details '%s'" % details)

    def write_audit_author(self):
        with self.loop("audit_author", ["name", "pdbx_ordinal"]) as l:
            l.write("'Pieper U' 1")
            l.write("'Webb B' 2")
            l.write("'Narayanan E' 3")
            l.write("'Sali A' 4")

    def write_citation(self):
        with self.loop("citation",
                ["id", "title", "journal_abbrev", "journal_volume",
                 "page_first", "page_last", "year", "pdbx_database_id_PubMed",
                 "pdbx_database_id_DOI"]) as l:
            # MODELLER paper
            l.write("1 'Comparative Protein Modelling by Satisfaction of "
                    "Spatial Restraints'\n'J Mol Biol' 234 779 815 1993 "
                    "8254673 10.1006/jmbi.1993.1626")
        with self.loop("citation_author",
                ["citation_id", "name", "ordinal"]) as l:
            l.write("1 'Sali A' 1")
            l.write("1 'Blundell T L' 2")

    def write_software(self, modpipe_version, modeller_version):
        with self.loop("software",
                ["pdbx_ordinal", "name", "classification",
                 "version", "type", "location", "citation_id"]) as l:
            l.write("1 ModPipe 'comparative modeling' %s "
                    "program https://salilab.org/modpipe/ ." % modpipe_version)
            l.write("2 MODELLER 'comparative modeling' %s program "
                    "https://salilab.org/modeller/ 1" % modeller_version)

        # Put each piece of software in its own group
        with self.loop("ma_software_group",
                ["ordinal_id", "group_id", "software_id"]) as l:
            for i in range(1,3):
                l.write("%d %d %d" % (i,i,i))

    def write_chem_comp(self):
        # Just assume all 20 standard amino acids are in the model
        with self.loop("chem_comp",
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

    def write_entity_details(self, sequence3):
        # entities for target (and template if we have alignment)
        with self.loop("entity",
                ["id", "type", "src_method", "pdbx_description"]) as l:
            if self.align:
                l.write("%d polymer man template" % self.template.entity_id)
            l.write("%d polymer man target" % self.target.entity_id)

        target_primary = "".join(three_to_one[x] for x in sequence3)

        with self.loop("entity_poly",
                ["entity_id", "type", "nstd_linkage",
                 "pdbx_seq_one_letter_code",
                 "pdbx_seq_one_letter_code_can"]) as l:
            if self.align:
                if self.align.target.primary != target_primary:
                    raise ValueError("Model sequence does not match target "
                            "sequence in alignment:",
                            target_primary, self.align.target.primary)
                p = self.align.template.primary
                l.write("%d polypeptide(L) no %s %s"
                        % (self.template.entity_id, p, p))
            l.write("%d polypeptide(L) no %s %s"
                    % (self.target.entity_id, target_primary, target_primary))

        with self.loop("entity_poly_seq",
                ["entity_id", "num", "mon_id", "hetero"]) as l:
            if self.align:
                p = self.align.template.primary
                for i, s in enumerate(p):
                    l.write("%d %d %s ."
                            % (self.template.entity_id, i+1, one_to_three[s]))
            for i, s in enumerate(sequence3):
                l.write("%d %d %s ." % (self.target.entity_id, i+1, s))

    def write_template_details(self, chain_id, pdb_beg, pdb_end, pdb_code):
        if not self.align:
            return
        # Define the identity transformation (id=1)
        self.print("#\n_ma_template_trans_matrix.id 1")
        for i in range(1,4):
            for j in range(1,4):
                self.print("_ma_template_trans_matrix.rot_matrix[%d][%d] %s"
                           % (j, i, "1.0" if i == j else "0.0"))
        for i in range(1,4):
            self.print("_ma_template_trans_matrix.tr_vector[%d] 0.0" % i)

        with self.loop("ma_template_details",
                  ["ordinal_id", "template_id", "template_origin",
                   "template_entity_type", "template_trans_matrix_id",
                   "template_data_id", "target_asym_id",
                   "template_label_asym_id", "template_label_entity_id",
                   "template_model_num"]) as l:
            # template structure is data_id=1
            # trans_matrix_id=1 is the identity transformation
            # model_num=1 because Modeller always uses the first PDB model
            l.write('1 1 "reference database" polymer 1 %d %s %s 1 1'
                    % (self.template.data_id, chain_id,
                       self.align.template.chain))

        with self.loop("ma_template_poly",
                  ["template_id", "seq_one_letter_code",
                   "seq_one_letter_code_can"]) as l:
            p = self.align.template.primary
            l.write("1 %s %s" % (p, p))

        if self.align:
            # template_id makes no sense if we have no alignment
            with self.loop("ma_template_poly_segment",
                      ["id", "template_id", "residue_number_begin",
                       "residue_number_end"]) as l:
                l.write("1 1 %d %d" % (pdb_beg, pdb_end))

        with self.loop("ma_template_ref_db_details",
                ["template_id", "db_name", "db_accession_code"]) as l:
            l.write("1 PDB %s" % pdb_code)

    def write_target_details(self, chain_id, sequence3):
        with self.loop("ma_target_entity",
                ["entity_id", "data_id", "origin"]) as l:
            l.write("%d %d ." % (self.target.entity_id, self.target.data_id))

        with self.loop("ma_target_entity_instance",
                ["asym_id", "entity_id", "details"]) as l:
            l.write("%s %d ." % (chain_id, self.target.entity_id))

        if self.align:
            # Cannot write a template segment ID without an alignment
            with self.loop("ma_target_template_poly_mapping",
                      ["id", "template_segment_id", "target_asym_id",
                       "target_seq_id_begin", "target_seq_id_end"]) as l:
                l.write("1 1 %s 1 %d" % (chain_id, len(sequence3)))

    def write_alignment(self, chain_id, evalue):
        if not self.align:
            return
        # Just one target-template alignment (one template, one chain) so this
        # table is pretty simple:
        with self.loop("ma_alignment_info",
                  ["alignment_id", "data_id", "software_id",
                   "alignment_length", "alignment_type",
                   "alignment_mode"]) as l:
            # ModPipe is software_id=1
            l.write('1 %d 1 %d "target-template pairwise alignment" global'
                    % (self.alignment.data_id, len(self.align.template.gapped)))

        with self.loop("ma_alignment_details",
                  ["ordinal_id", "alignment_id", "template_segment_id",
                   "target_asym_id", "score_type", "score_value"]) as l:
            l.write("1 1 1 %s 'BLAST e-value' %s" % (chain_id, evalue))

        with self.loop("ma_alignment",
                  ["ordinal_id", "alignment_id", "target_template_flag",
                   "sequence"]) as l:
            # Template (flag=2)
            l.write("1 1 2 %s" % self.align.template.gapped)
            # Target (flag=1)
            l.write("2 1 1 %s" % self.align.target.gapped)

    def write_assembly(self, chain_id, sequence3):
        with self.loop('ma_struct_assembly',
                ['ordinal_id', 'assembly_id', 'entity_id', 'asym_id',
                 'seq_id_begin', 'seq_id_end']) as l:
            # Simple assembly of a single chain
            l.write("1 1 %d %s 1 %d"
                    % (self.target.entity_id, chain_id, len(sequence3)))

    def write_data(self):
        with self.loop("ma_data", ["id", "name", "content_type"]) as l:
            l.write("%d 'Template Structure' 'template structure'"
                    % self.template.data_id)
            l.write("%d 'Target Sequence' target" % self.target.data_id)
            l.write("%d 'Target Template Alignment' "
                    "'target-template alignment'" % self.alignment.data_id)
            l.write("%d 'Target Structure' 'model coordinates'"
                    % self.coord.data_id)

        # Put each data item in its own group
        with self.loop("ma_data_group",
                ["ordinal_id", "group_id", "data_id"]) as l:
            for i in range(1, 5):
                l.write("%d %d %d" % (i,i,i))

    def write_protocol(self):
        with self.loop('ma_protocol_step',
                ['ordinal_id', 'protocol_id', 'step_id', 'method_type',
                 'step_name', 'software_group_id', 'input_data_group_id',
                 'output_data_group_id']) as l:
            # step 1, template search
            # template source is the ModPipe fold assignment method
            # ModPipe is software_id=1
            # takes as input the template
            # makes as output the alignment
            l.write("1 1 1 'template search' 'ModPipe %s' 1 %d %d"
                    % (self.align.template.method if self.align else '.',
                       self.template.data_id, self.alignment.data_id))

            # step 2, modeling
            # MODELLER is software_id=2
            # takes as input the alignment
            # makes as output the coordinates
            l.write("2 1 2 'modeling' . 2 %d %d"
                    % (self.alignment.data_id, self.coord.data_id))

            # step 3, model selection
            # takes as input the coordinates and returns them unchanged
            l.write("3 1 3 'model selection' . 1 %d %d"
                    % (self.coord.data_id, self.coord.data_id))

    def write_scores(self, tsvmod_method, tsvmod_rmsd, tsvmod_no35, ga341,
                     zdope, mpqs):
        if not mpqs:
            return
        with self.loop('ma_qa_metric',
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

        with self.loop('ma_qa_metric_global',
                ['ordinal_id', 'model_id', 'metric_id', 'value']) as l:
            l.write("1 1 1 %s" % mpqs)
            l.write("2 1 2 %s" % zdope)
            if tsvmod_rmsd:
                l.write("3 1 3 %s" % tsvmod_rmsd)
                l.write("4 1 4 %s" % tsvmod_no35)

    def write_model_list(self):
        with self.loop('ma_model_list',
                ['ordinal_id', 'model_id', 'model_group_id', 'model_name',
                 'model_group_name', 'assembly_id', 'data_id',
                 'model_type']) as l:
            l.write("1 1 1 'Selected model' . 1 %d 'Homology model'"
                    % self.coord.data_id)

    def write_asym(self, chain_id):
        with self.loop('struct_asym', ['id', 'entity_id', 'details']) as l:
            l.write("%s %d ?" % (chain_id, self.target.entity_id))

    def write_atom_site(self, chain_id, atoms, resnum_begin, resnum_end):
        elements = set()
        auth_seqid = resnum_begin
        entity_id = self.target.entity_id
        seqid = 1
        ordinal = 1
        pdb_resnum = None
        with self.loop('atom_site',
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
                l.write("%s %s %s %s %s %d %d %s %s . %s %s %s %s %s %d %d"
                        % (group_pdb, element, atmnam, resnam, chain_id, seqid,
                           auth_seqid, inscode, chain_id, x, y, z, occ, tfac,
                           entity_id, ordinal))
                ordinal += 1
        assert auth_seqid == resnum_end

        with CifLoop(fh, 'atom_type', ['symbol']) as l:
            for element in sorted(elements):
                l.write(element)


class Structure:
    """Handle read of PDB structure and write of mmCIF"""

    def _read_pdb(self, fh):
        self.remarks = {}
        self.expdta = None
        self.title = None
        self.modpipe_version = None
        self.atoms = []

        for line in fh:
            # Handle standard ModBase headers
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
        """Get PDB sequence as a sequence of 3-letter residue names"""
        resnum = None
        for a in self.atoms:
            this_resnum = a[22:26]  # residue sequence number
            if this_resnum != resnum:
                yield a[17:20].strip()  # residue name
                resnum = this_resnum

    def write_mmcif(self, fh, align):
        """Write current structure out to a mmCIF file handle"""
        # mmCIF models must always have a chain ID; older ModBase PDB models had
        # a blank ID
        chain_id = self.chain_id or 'A'
        sequence3 = list(self.get_sequence3())
        modeller_version = self.get_modeller_version() or '?'
        if align:
            align = Alignment(align)

        c = CifWriter(fh, align)
        c.write_header(self.remarks['MODPIPE MODEL ID'], self.title)
        c.write_exptl(self.remarks['MODPIPE MODEL ID'], self.expdta)
        c.write_audit_author()
        c.write_citation()
        c.write_software(self.modpipe_version, modeller_version)
        c.write_chem_comp()
        c.write_entity_details(sequence3)
        c.write_template_details(chain_id, int(self.remarks['TEMPLATE BEGIN']),
                int(self.remarks['TEMPLATE END']), self.remarks['TEMPLATE PDB'])
        c.write_target_details(chain_id, sequence3)
        c.write_alignment(chain_id, self.remarks['EVALUE'])
        c.write_assembly(chain_id, sequence3)
        c.write_data()
        c.write_protocol()
        c.write_scores(self.remarks.get('TSVMOD METHOD'),
                self.remarks.get('TSVMOD RMSD'),
                self.remarks.get('TSVMOD NO35'),
                self.remarks.get('GA341 SCORE'),
                self.remarks.get('zDOPE SCORE'),
                self.remarks.get('MPQS'))
        c.write_model_list()
        c.write_asym(chain_id)
        c.write_atom_site(chain_id, self.atoms,
                int(self.remarks['TARGET BEGIN']),
                int(self.remarks['TARGET END']))


def read_pdb(fh):
    """Read PDB file from filehandle and return a new Structure"""
    s = Structure()
    s._read_pdb(fh)
    return s


if __name__ == '__main__':
    import argparse
    a = argparse.ArgumentParser(
            description="Utility to convert ModBase PDB files to mmCIF")
    a.add_argument("-a", "--align", metavar="FILE", help="Input alignment file")
    a.add_argument("pdb", help="Input PDB file")
    a.add_argument("mmcif", help="Output mmCIF file")
    args = a.parse_args()

    with open(args.pdb) as fh:
        s = read_pdb(fh)
    with open(args.mmcif, 'w') as fh:
        s.write_mmcif(fh, args.align)
