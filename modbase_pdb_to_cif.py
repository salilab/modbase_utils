#!/usr/bin/python3

import re
import os
import collections
import operator
import gzip
import itertools
import ihm.reader
import ihm.format
import ihm.format_bcif

# Single sequence in a Modeller alignment
Sequence = collections.namedtuple(
    "Sequence", ["seqtyp", "chain", "method", "gapped", "primary"])


# Reference sequence database
SequenceDB = collections.namedtuple(
    "SequenceDB", ["name", "code", "accession"])


# Mapping between one-letter codes and PDB names
three_to_one = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

one_to_three = {val: key for key, val in three_to_one.items()}


def split_resnum(resnum):
    """Split a residue number into number and insertion code (or None)"""
    m = re.match(r'([\d-]+)(.*)$', resnum)
    return m.group(1), m.group(2) or None


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
        return Sequence(
            seqtyp=header[0], chain=header[3], method=header[7],
            gapped=gapped, primary=gapped.replace('-', ''))


class CifWriter:
    def __init__(self, writer, align):
        class CifID:
            pass

        self.writer = writer
        self.align = align
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

    def flush(self):
        self.writer.flush()

    def loop(self, category, keys):
        return self.writer.loop(category, keys)

    def write_header(self, model_id, title):
        self.writer.start_block("model_%s" % model_id)
        entry_id = "modbase_model_%s" % model_id
        with self.writer.category("_entry") as lp:
            lp.write(id=entry_id)
        with self.writer.category("_struct") as lp:
            lp.write(entry_id=entry_id, title=title)

    def write_exptl(self, model_id, expdta):
        if expdta.startswith('THEORETICAL MODEL, '):
            details = expdta[19:]
            with self.writer.category("_exptl") as cat:
                cat.write(entry_id="modbase_model_%s" % model_id,
                          method='THEORETICAL MODEL', details=details)

    def write_audit_author(self):
        ordinal = itertools.count(1)
        with self.loop("_audit_author", ["name", "pdbx_ordinal"]) as lp:
            lp.write(name='Pieper U', pdbx_ordinal=next(ordinal))
            lp.write(name='Webb B', pdbx_ordinal=next(ordinal))
            lp.write(name='Narayanan E', pdbx_ordinal=next(ordinal))
            lp.write(name='Sali A', pdbx_ordinal=next(ordinal))

    def write_citation(self):
        with self.loop(
                "_citation",
                ["id", "title", "journal_abbrev", "journal_volume",
                 "page_first", "page_last", "year", "pdbx_database_id_PubMed",
                 "pdbx_database_id_DOI"]) as lp:
            # MODELLER paper
            lp.write(id=1,
                     title="Comparative Protein Modelling by Satisfaction of "
                           "Spatial Restraints",
                     journal_abbrev='J Mol Biol', journal_volume=234,
                     page_first=779, page_last=815, year=1993,
                     pdbx_database_id_PubMed=8254673,
                     pdbx_database_id_DOI="10.1006/jmbi.1993.1626")
        ordinal = itertools.count(1)
        with self.loop(
                "_citation_author", ["citation_id", "name", "ordinal"]) as lp:
            lp.write(citation_id=1, name='Sali A', ordinal=next(ordinal))
            lp.write(citation_id=1, name='Blundell T L', ordinal=next(ordinal))

    def write_software(self, modpipe_version, modeller_version):
        ordinal = itertools.count(1)
        with self.loop(
                "_software",
                ["pdbx_ordinal", "name", "classification",
                 "version", "type", "location", "citation_id"]) as lp:
            lp.write(pdbx_ordinal=next(ordinal),
                     name="ModPipe", classification='comparative modeling',
                     version=modpipe_version, type="program",
                     location="https://salilab.org/modpipe/",
                     citation_id=None)
            lp.write(pdbx_ordinal=next(ordinal),
                     name="MODELLER", classification='comparative modeling',
                     version=modeller_version, type="program",
                     location="https://salilab.org/modeller/",
                     citation_id=1)

        # Put each piece of software in its own group
        with self.loop(
                "_ma_software_group",
                ["ordinal_id", "group_id", "software_id"]) as lp:
            for i in range(1, 3):
                lp.write(ordinal_id=i, group_id=i, software_id=i)

    def write_chem_comp(self, sequence3):
        seen_one_letter = set(three_to_one[x] for x in sequence3)
        if self.align:
            seen_one_letter.update(self.align.template.primary)
        al = ihm.LPeptideAlphabet()
        comps = [al[code] for code in seen_one_letter]

        with self.loop(
                "_chem_comp",
                ["id", "type", "name", "formula", "formula_weight"]) as lp:
            for cc in sorted(comps, key=operator.attrgetter('id')):
                lp.write(id=cc.id, type=cc.type, name=cc.name,
                         formula=cc.formula, formula_weight=cc.formula_weight)

    def write_entity_details(self, sequence3):
        # entities for target (and template if we have alignment)
        with self.loop(
                "_entity",
                ["id", "type", "src_method", "pdbx_description"]) as lp:
            if self.align:
                lp.write(id=self.template.entity_id, type="polymer",
                         src_method="man", pdbx_description="template")
            lp.write(id=self.target.entity_id, type="polymer",
                     src_method="man", pdbx_description="target")

        target_primary = "".join(three_to_one[x] for x in sequence3)

        with self.loop(
                "_entity_poly",
                ["entity_id", "type", "nstd_linkage",
                 "pdbx_seq_one_letter_code",
                 "pdbx_seq_one_letter_code_can"]) as lp:
            if self.align:
                if self.align.target.primary != target_primary:
                    raise ValueError(
                        "Model sequence does not match target "
                        "sequence in alignment:",
                        target_primary, self.align.target.primary)
                p = self.align.template.primary
                lp.write(entity_id=self.template.entity_id,
                         type="polypeptide(L)", nstd_linkage="no",
                         pdbx_seq_one_letter_code=p,
                         pdbx_seq_one_letter_code_can=p)
            lp.write(entity_id=self.target.entity_id,
                     type="polypeptide(L)", nstd_linkage="no",
                     pdbx_seq_one_letter_code=target_primary,
                     pdbx_seq_one_letter_code_can=target_primary)

        with self.loop(
                "_entity_poly_seq",
                ["entity_id", "num", "mon_id", "hetero"]) as lp:
            if self.align:
                p = self.align.template.primary
                for i, s in enumerate(p):
                    lp.write(entity_id=self.template.entity_id, num=i+1,
                             mon_id=one_to_three[s], hetero=None)
            for i, s in enumerate(sequence3):
                lp.write(entity_id=self.target.entity_id, num=i+1,
                         mon_id=s, hetero=None)

    def write_template_details(self, chain_id, tmpbeg, tmpend, tmpasym,
                               pdb_code):
        if not self.align:
            return
        # Define the identity transformation (id=1)
        with self.writer.loop(
                "_ma_template_trans_matrix",
                ["id",
                 "rot_matrix[1][1]", "rot_matrix[2][1]", "rot_matrix[3][1]",
                 "rot_matrix[1][2]", "rot_matrix[2][2]", "rot_matrix[3][2]",
                 "rot_matrix[1][3]", "rot_matrix[2][3]", "rot_matrix[3][3]",
                 "tr_vector[1]", "tr_vector[2]", "tr_vector[3]"]) as lp:
            lp.write(id=1, rot_matrix11="1.0", rot_matrix22="1.0",
                     rot_matrix33="1.0", rot_matrix12="0.0",
                     rot_matrix13="0.0", rot_matrix21="0.0",
                     rot_matrix23="0.0", rot_matrix31="0.0",
                     rot_matrix32="0.0", tr_vector1="0.0",
                     tr_vector2="0.0", tr_vector3="0.0")

        with self.loop(
                "_ma_template_details",
                ["ordinal_id", "template_id", "template_origin",
                 "template_entity_type", "template_trans_matrix_id",
                 "template_data_id", "target_asym_id",
                 "template_label_asym_id", "template_label_entity_id",
                 "template_model_num"]) as lp:
            # template structure is data_id=1
            # trans_matrix_id=1 is the identity transformation
            # model_num=1 because Modeller always uses the first PDB model
            lp.write(ordinal_id=1, template_id=1,
                     template_origin="reference database",
                     template_entity_type="polymer",
                     template_trans_matrix_id=1,
                     template_data_id=self.template.data_id,
                     target_asym_id=chain_id, template_label_asym_id=tmpasym,
                     template_label_entity_id=1, template_model_num=1)

        with self.loop(
                "_ma_template_poly",
                ["template_id", "seq_one_letter_code",
                 "seq_one_letter_code_can"]) as lp:
            p = self.align.template.primary
            lp.write(template_id=1, seq_one_letter_code=p,
                     seq_one_letter_code_can=p)

        if self.align:
            # template_id makes no sense if we have no alignment
            with self.loop(
                    "_ma_template_poly_segment",
                    ["id", "template_id", "residue_number_begin",
                     "residue_number_end"]) as lp:
                lp.write(id=1, template_id=1, residue_number_begin=tmpbeg,
                         residue_number_end=tmpend)

        with self.loop(
                "_ma_template_ref_db_details",
                ["template_id", "db_name", "db_accession_code"]) as lp:
            lp.write(template_id=1, db_name="PDB", db_accession_code=pdb_code)

    def write_target_details(self, chain_id, sequence3, seqdb, tgtbeg, tgtend):
        with self.loop(
                "_ma_target_entity", ["entity_id", "data_id", "origin"]) as lp:
            lp.write(entity_id=self.target.entity_id,
                     data_id=self.target.data_id, origin=None)

        with self.loop(
                "_ma_target_entity_instance",
                ["asym_id", "entity_id", "details"]) as lp:
            lp.write(asym_id=chain_id, entity_id=self.target.entity_id,
                     details=None)

        if self.align:
            # Cannot write a template segment ID without an alignment
            with self.loop(
                    "_ma_target_template_poly_mapping",
                    ["id", "template_segment_id", "target_asym_id",
                     "target_seq_id_begin", "target_seq_id_end"]) as lp:
                lp.write(id=1, template_segment_id=1, target_asym_id=chain_id,
                         target_seq_id_begin=1,
                         target_seq_id_end=len(sequence3))

        with self.loop(
                "_ma_target_ref_db_details",
                ["target_entity_id", "db_name", "db_name_other_details",
                 "db_code", "db_accession", "seq_db_isoform",
                 "seq_db_align_begin", "seq_db_align_end"]) as lp:
            for db in seqdb:
                if db.name == 'UniProt':
                    lp.write(target_entity_id=self.target.entity_id,
                             db_name="UNP", db_name_other_details=None,
                             db_code=db.code, db_accession=db.accession,
                             seq_db_isoform=ihm.unknown,
                             seq_db_align_begin=tgtbeg,
                             seq_db_align_end=tgtend)
                elif db.name in ('RefSeq', 'PlasmoDB'):
                    lp.write(target_entity_id=self.target.entity_id,
                             db_name="Other", db_name_other_details=db.name,
                             db_code=db.code, db_accession=db.accession,
                             seq_db_isoform=ihm.unknown,
                             seq_db_align_begin=tgtbeg,
                             seq_db_align_end=tgtend)

    def write_alignment(self, chain_id, evalue):
        if not self.align:
            return
        # Just one target-template alignment (one template, one chain) so this
        # table is pretty simple:
        with self.loop(
                "_ma_alignment_info",
                ["alignment_id", "data_id", "software_group_id",
                 "alignment_length", "alignment_type",
                 "alignment_mode"]) as lp:
            # ModPipe is software_group_id=1
            lp.write(alignment_id=1, data_id=self.alignment.data_id,
                     software_group_id=1,
                     alignment_length=len(self.align.template.gapped),
                     alignment_type="target-template pairwise alignment",
                     alignment_mode="global")

        with self.loop(
                "_ma_alignment_details",
                ["ordinal_id", "alignment_id", "template_segment_id",
                 "target_asym_id", "score_type", "score_value"]) as lp:
            lp.write(ordinal_id=1, alignment_id=1, template_segment_id=1,
                     target_asym_id=chain_id, score_type='BLAST e-value',
                     score_value=evalue)

        with self.loop(
                "_ma_alignment",
                ["ordinal_id", "alignment_id", "target_template_flag",
                 "sequence"]) as lp:
            ordinal = itertools.count(1)
            # Template (flag=2)
            lp.write(ordinal_id=next(ordinal), alignment_id=1,
                     target_template_flag=2,
                     sequence=self.align.template.gapped)
            # Target (flag=1)
            lp.write(ordinal_id=next(ordinal), alignment_id=1,
                     target_template_flag=1,
                     sequence=self.align.target.gapped)

    def write_assembly(self, chain_id, sequence3):
        with self.loop(
                '_ma_struct_assembly',
                ['ordinal_id', 'assembly_id', 'entity_id', 'asym_id',
                 'seq_id_begin', 'seq_id_end']) as lp:
            # Simple assembly of a single chain
            lp.write(ordinal_id=1, assembly_id=1,
                     entity_id=self.target.entity_id, asym_id=chain_id,
                     seq_id_begin=1, seq_id_end=len(sequence3))

    def write_data(self):
        with self.loop("_ma_data", ["id", "name", "content_type"]) as lp:
            lp.write(id=self.template.data_id, name='Template Structure',
                     content_type='template structure')
            lp.write(id=self.target.data_id, name='Target Sequence',
                     content_type='target')
            lp.write(id=self.alignment.data_id,
                     name='Target Template Alignment',
                     content_type='target-template alignment')
            lp.write(id=self.coord.data_id,
                     name='Target Structure',
                     content_type='model coordinates')

        # Put each data item in its own group
        with self.loop(
                "_ma_data_group", ["ordinal_id", "group_id", "data_id"]) as lp:
            for i in range(1, 5):
                lp.write(ordinal_id=i, group_id=i, data_id=i)

    def write_protocol(self):
        with self.loop(
                '_ma_protocol_step',
                ['ordinal_id', 'protocol_id', 'step_id', 'method_type',
                 'step_name', 'software_group_id', 'input_data_group_id',
                 'output_data_group_id']) as lp:
            ordinal = itertools.count(1)
            # step 1, template search
            # template source is the ModPipe fold assignment method
            # ModPipe is software_id=1
            # takes as input the template
            # makes as output the alignment
            mth = " %s" % self.align.template.method if self.align else ''
            lp.write(ordinal_id=next(ordinal), protocol_id=1, step_id=1,
                     method_type='template search',
                     step_name='ModPipe%s' % mth, software_group_id=1,
                     input_data_group_id=self.template.data_id,
                     output_data_group_id=self.alignment.data_id)

            # step 2, modeling
            # MODELLER is software_id=2
            # takes as input the alignment
            # makes as output the coordinates
            lp.write(ordinal_id=next(ordinal), protocol_id=1, step_id=2,
                     method_type='modeling',
                     step_name=None, software_group_id=2,
                     input_data_group_id=self.alignment.data_id,
                     output_data_group_id=self.coord.data_id)

            # step 3, model selection
            # takes as input the coordinates and returns them unchanged
            lp.write(ordinal_id=next(ordinal), protocol_id=1, step_id=3,
                     method_type='model selection',
                     step_name=None, software_group_id=1,
                     input_data_group_id=self.coord.data_id,
                     output_data_group_id=self.coord.data_id)

    def write_scores(self, tsvmod_method, tsvmod_rmsd, tsvmod_no35, ga341,
                     zdope, mpqs):
        if not mpqs:
            return
        ordinal = itertools.count(1)
        with self.loop(
                '_ma_qa_metric',
                ['id', 'name', 'description', 'type', 'mode',
                 'type_other_details', 'software_group_id']) as lp:
            # ModPipe is software_id=1
            lp.write(id=next(ordinal), name='MPQS',
                     description='ModPipe Quality Score',
                     type="other", mode="global",
                     type_other_details="composite score, values >1.1 are "
                                        "considered reliable",
                     software_group_id=1)

            # MODELLER is software_id=2
            lp.write(id=next(ordinal), name="zDOPE",
                     description='Normalized DOPE',
                     type="zscore", mode="global", software_group_id=2)

            if tsvmod_rmsd:
                lp.write(id=next(ordinal), name='TSVMod RMSD',
                         description='TSVMod predicted RMSD (%s)'
                                     % tsvmod_method,
                         type="distance", mode="global")
                lp.write(id=next(ordinal), name='TSVMod NO35',
                         description="TSVMod predicted native "
                                     "overlap (%s)" % tsvmod_method,
                         type="other", mode="global")

        with self.loop(
                '_ma_qa_metric_global',
                ['ordinal_id', 'model_id', 'metric_id', 'metric_value']) as lp:
            ordinal = itertools.count(1)
            lp.write(ordinal_id=next(ordinal), model_id=1, metric_id=1,
                     metric_value=mpqs)
            lp.write(ordinal_id=next(ordinal), model_id=1, metric_id=2,
                     metric_value=zdope)
            if tsvmod_rmsd:
                lp.write(ordinal_id=next(ordinal), model_id=1, metric_id=3,
                         metric_value=tsvmod_rmsd)
                lp.write(ordinal_id=next(ordinal), model_id=1, metric_id=4,
                         metric_value=tsvmod_no35)

    def write_model_list(self):
        with self.loop(
                '_ma_model_list',
                ['ordinal_id', 'model_id', 'model_group_id', 'model_name',
                 'model_group_name', 'assembly_id', 'data_id',
                 'model_type']) as lp:
            lp.write(ordinal_id=1, model_id=1, model_group_id=1,
                     model_name='Selected model', assembly_id=1,
                     data_id=self.coord.data_id, model_type='Homology model')

    def write_asym(self, chain_id):
        with self.loop('_struct_asym', ['id', 'entity_id', 'details']) as lp:
            lp.write(id=chain_id, entity_id=self.target.entity_id,
                     details=ihm.unknown)

    def write_seq_scheme(self, chain_id, sequence3, tgtbeg, tgtend):
        assert len(sequence3) == tgtend - tgtbeg + 1
        entity_id = self.target.entity_id
        with self.loop(
            '_pdbx_poly_seq_scheme',
            ['asym_id', 'entity_id', 'seq_id', 'mon_id', 'pdb_seq_num',
             'auth_seq_num', 'pdb_mon_id', 'auth_mon_id',
             'pdb_strand_id']) as lp:
            for i, s in enumerate(sequence3):
                seqid = 1 + i
                auth_seqid = tgtbeg + i
                lp.write(asym_id=chain_id, entity_id=entity_id, seq_id=seqid,
                         mon_id=s, pdb_seq_num=auth_seqid,
                         auth_seq_num=auth_seqid, pdb_mon_id=s, auth_mon_id=s,
                         pdb_strand_id=chain_id)

    def write_atom_site(self, chain_id, atoms, resnum_begin, resnum_end):
        elements = set()
        auth_seqid = resnum_begin
        entity_id = self.target.entity_id
        seqid = 1
        ordinal = itertools.count(1)
        pdb_resnum = None
        with self.loop(
            '_atom_site',
            ['group_PDB', 'type_symbol', 'label_atom_id',
             'label_comp_id', 'label_asym_id', 'label_seq_id',
             'auth_seq_id', 'pdbx_PDB_ins_code', 'auth_asym_id',
             'label_alt_id', 'Cartn_x', 'Cartn_y', 'Cartn_z',
             'occupancy', 'B_iso_or_equiv', 'label_entity_id',
             'pdbx_PDB_model_num', 'id']) as lp:
            for a in atoms:
                # Detect new residue if PDB resnum changed
                pdb_this_resnum = a[22:26]
                if pdb_resnum is not None and pdb_this_resnum != pdb_resnum:
                    auth_seqid += 1
                    seqid += 1
                pdb_resnum = pdb_this_resnum
                inscode = a[26:27].strip() or ihm.unknown
                group_pdb = a[:6].strip()
                element = a[76:78].strip() or ihm.unknown
                elements.add(element)
                atmnam = a[12:16].strip()
                resnam = a[17:20].strip()
                x = a[30:38].strip()
                y = a[38:46].strip()
                z = a[46:54].strip()
                occ = a[54:60].strip()
                tfac = a[60:66].strip()
                lp.write(group_PDB=group_pdb, type_symbol=element,
                         label_atom_id=atmnam, label_comp_id=resnam,
                         label_asym_id=chain_id, label_seq_id=seqid,
                         auth_seq_id=auth_seqid, pdbx_PDB_ins_code=inscode,
                         auth_asym_id=chain_id, label_alt_id=None, Cartn_x=x,
                         Cartn_y=y, Cartn_z=z, occupancy=occ,
                         B_iso_or_equiv=tfac, label_entity_id=entity_id,
                         pdbx_PDB_model_num=1,
                         id=next(ordinal))
        assert auth_seqid == resnum_end

        with self.loop('_atom_type', ['symbol']) as lp:
            for element in sorted(elements):
                lp.write(symbol=element)


class _PolySeqSchemeHandler(ihm.reader.Handler):
    """Read pdbx_poly_seq_scheme table and map PDB to mmCIF numbering"""
    def __init__(self, m):
        self.m = m

    def __call__(self, asym_id, seq_id, auth_seq_num, pdb_ins_code,
                 pdb_strand_id):
        mk = (pdb_strand_id, auth_seq_num, pdb_ins_code)
        if mk in self.m:
            self.m[mk] = (asym_id, seq_id)


class Repository:
    """Point to a directory containing a mirror of PDB in mmCIF format"""
    def __init__(self, topdir):
        self.topdir = topdir

    def open_mmcif(self, pdb_code):
        """Given a PDB code, return a file handle to the corresponding mmCIF"""
        code = pdb_code.lower()
        fname = os.path.join(self.topdir, code[1:3], code + '.cif.gz')
        return gzip.open(fname, 'rt', encoding='latin1')

    def map_ranges(self, fh, ranges):
        """Map a list of PDB (chain, resnum_start, resnum_end) tuples to the
           corresponding mmCIF (asym_id, seq_id_start, seq_id_end) values
           and return. This is done by reading the pdbx_poly_seq_scheme
           table in the given mmCIF file."""
        m = collections.OrderedDict()
        for r in ranges:
            # If we can't find the PDB numbering, return it unchanged
            m[r] = (r[0], r[1])
        h = _PolySeqSchemeHandler(m)
        r = ihm.format.CifReader(fh, {'_pdbx_poly_seq_scheme': h})
        r.read_file()
        return m.values()


class Structure:
    """Handle read of PDB structure and write of mmCIF"""

    def __init__(self, repo):
        self.repo = repo

    def _read_pdb(self, fh):
        self.remarks = {}
        self.seqdb = []
        self.expdta = None
        self.title = None
        self.modpipe_version = None
        self.atoms = []

        for line in fh:
            # Handle standard ModBase headers
            if line.startswith('REMARK 220 SEQDB:'):
                val = [x.strip() for x in line[17:].split()]
                if len(val) == 3:
                    self.seqdb.append(SequenceDB(
                        name=val[0], accession=val[1], code=val[2]))
                elif len(val) == 2 and val[0] == 'RefSeq' and '.' in val[1]:
                    self.seqdb.append(SequenceDB(
                        name=val[0], accession=val[1].split('.', 1)[0],
                        code=val[1]))
                elif len(val) == 2 and val[0] in ('PlasmoDB',):
                    self.seqdb.append(SequenceDB(
                        name=val[0], accession=val[1], code='.'))
            elif line.startswith('REMARK') and line.count(':') == 1:
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
            m = re.search(r'MODELLER\s+(\S+)', self.expdta)
            if m:
                return m.group(1)

    def get_mmcif_template_info(self, pdb_beg, pdb_end, pdb_chain, pdb_code):
        """Given PDB ("author-provided") template information, map to
           mmCIF seq_id range and asym_id, and return."""
        # Split TEMPLATE BEGIN/END records into residue number and
        # insertion code
        pdb_beg, pdb_beg_ins = split_resnum(pdb_beg)
        pdb_end, pdb_end_ins = split_resnum(pdb_end)
        pdb_ranges = [(pdb_chain, pdb_beg, pdb_beg_ins),
                      (pdb_chain, pdb_end, pdb_end_ins)]
        # Open the mmCIF file and map PDB ranges to mmCIF
        with self.repo.open_mmcif(pdb_code) as fh:
            cif_ranges = list(self.repo.map_ranges(fh, pdb_ranges))
        # asym_id of start and end should be the same
        assert(cif_ranges[0][0] == cif_ranges[1][0])
        return(int(cif_ranges[0][1]), int(cif_ranges[1][1]), cif_ranges[0][0])

    def get_sequence3(self):
        """Get PDB sequence as a sequence of 3-letter residue names"""
        resnum = None
        for a in self.atoms:
            this_resnum = a[22:26]  # residue sequence number
            if this_resnum != resnum:
                yield a[17:20].strip()  # residue name
                resnum = this_resnum

    def write_mmcif(self, writer, align):
        """Write current structure out to a mmCIF file handle"""
        # mmCIF models must always have a chain ID; older ModBase PDB models
        # had a blank ID
        chain_id = self.chain_id.strip() or 'A'
        sequence3 = list(self.get_sequence3())
        modeller_version = self.get_modeller_version() or '?'
        if align:
            align = Alignment(align)

        c = CifWriter(writer, align)
        c.write_header(self.remarks['MODPIPE MODEL ID'], self.title)
        c.write_exptl(self.remarks['MODPIPE MODEL ID'], self.expdta)
        c.write_audit_author()
        c.write_citation()
        c.write_software(self.modpipe_version, modeller_version)
        c.write_chem_comp(sequence3)
        c.write_entity_details(sequence3)
        template_pdb = self.remarks['TEMPLATE PDB']
        tmpbeg, tmpend, tmpasym = self.get_mmcif_template_info(
            self.remarks['TEMPLATE BEGIN'], self.remarks['TEMPLATE END'],
            self.remarks['TEMPLATE CHAIN'], template_pdb)
        c.write_template_details(
            chain_id, tmpbeg, tmpend, tmpasym, template_pdb)
        tgtbeg = int(self.remarks['TARGET BEGIN'])
        tgtend = int(self.remarks['TARGET END'])
        c.write_target_details(chain_id, sequence3, self.seqdb, tgtbeg, tgtend)
        c.write_alignment(chain_id, self.remarks['EVALUE'])
        c.write_assembly(chain_id, sequence3)
        c.write_data()
        c.write_protocol()
        c.write_scores(
            self.remarks.get('TSVMOD METHOD'), self.remarks.get('TSVMOD RMSD'),
            self.remarks.get('TSVMOD NO35'),
            self.remarks.get('GA341 SCORE', self.remarks.get('MODEL SCORE')),
            self.remarks.get('zDOPE SCORE', self.remarks.get('ZDOPE SCORE')),
            self.remarks.get('MPQS',
                             self.remarks.get('MODPIPE QUALITY SCORE')))
        c.write_model_list()
        c.write_asym(chain_id)
        c.write_seq_scheme(chain_id, sequence3, tgtbeg, tgtend)
        c.write_atom_site(chain_id, self.atoms, tgtbeg, tgtend)
        c.flush()


def read_pdb(fh, repo):
    """Read PDB file from filehandle and return a new Structure"""
    s = Structure(repo)
    s._read_pdb(fh)
    return s


if __name__ == '__main__':
    import argparse
    a = argparse.ArgumentParser(
            description="Utility to convert ModBase PDB files to mmCIF",
            epilog="""
Convert a PDB file, downloaded from ModBase, to mmCIF format. This should
preserve all information in the PDB file. If the corresponding alignment is
also provided (-a flag) alignment information is added to the mmCIF (in which
case the resulting mmCIF file should match that downloaded directly from
ModBase).
""")
    a.add_argument("-a", "--align", metavar="FILE",
                   help="Input alignment file")
    a.add_argument("-r", "--repo", default='.',
                   help="Directory containing repository of mmCIF files")
    a.add_argument("pdb", help="Input PDB file")
    a.add_argument("mmcif", help="Output mmCIF file")
    args = a.parse_args()

    r = Repository(args.repo)
    with open(args.pdb) as fh:
        s = read_pdb(fh, r)
    if args.mmcif.endswith('.bcif'):
        mode, writer = "wb", ihm.format_bcif.BinaryCifWriter
    else:
        mode, writer = "w", ihm.format.CifWriter
    with open(args.mmcif, mode) as fh:
        s.write_mmcif(writer(fh), args.align)
