#!/usr/bin/python3

import re
import os
import sys
import collections
import operator
import gzip
import itertools
import ihm.reader
import ihm.citations
import ma
import ma.dumper
import ma.model
import ma.alignment
import ma.reference
from ma.alignment import ShorterSequenceIdentity as SequenceIdentity
import ma.protocol

class RefSeq(ma.reference.TargetReference):
    other_details = "RefSeq"


class PlasmoDB(ma.reference.TargetReference):
    other_details = "PlasmoDB"


class MPQSMetricType(ma.qa_metric.MetricType):
    other_details = "composite score, values >1.1 are considered reliable"


# Single sequence in a Modeller alignment
Sequence = collections.namedtuple(
    "Sequence", ["seqtyp", "chain", "method", "gapped", "primary",
                 "primary_can"])


# Reference sequence database
SequenceDB = collections.namedtuple(
    "SequenceDB", ["name", "code", "accession"])


# Mapping between one-letter codes and PDB names
three_to_one = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
        'GLX': 'Z', 'ASX': 'B'
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
        # "Canonical" primary sequence is always a sequence of one-letter
        # codes; regular primary sequence is 1-letter for standard amino
        # acids, but can be longer for any non-standard residues (currently
        # only UNK is handled here, assuming X always means UNK in the
        # template).
        primary_can = gapped.replace('-', '')
        primary = ['UNK' if x == 'X' else x for x in primary_can]
        return Sequence(
            seqtyp=header[0], chain=header[3], method=header[7],
            gapped=gapped, primary=primary, primary_can=primary_can)


class _PolySeqSchemeHandler(ihm.reader.Handler):
    """Read pdbx_poly_seq_scheme table and map PDB to mmCIF numbering"""
    def __init__(self, m):
        self.m = m

    def __call__(self, asym_id, seq_id, pdb_seq_num, pdb_ins_code,
                 pdb_strand_id):
        mk = (pdb_strand_id, pdb_seq_num, pdb_ins_code)
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
            elif line.startswith('REMARK   6 GENERATED BY MODPIPE VERSION '):
                self.modpipe_version = line[40:].strip()
            elif line.startswith('REMARK') and line.count(':') == 1:
                key, val = [x.strip() for x in line[11:].split(':')]
                self.remarks[key] = val
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
            m = re.search(r'MODELLER\s+(Version\s+)?(\S+)', self.expdta)
            if m:
                return m.group(2)

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
        # Handle start==end
        if len(cif_ranges) == 1:
            cif_ranges = [cif_ranges[0], cif_ranges[0]]
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

    def get_system(self, align):
        """Create and return an ma.System object"""
        if align:
            align = Alignment(align)
        tgt_primary = "".join(three_to_one[x] for x in self.get_sequence3())

        s = ma.System(title=self.title,
                      id='model_' + self.remarks['MODPIPE MODEL ID'])
        s.authors.extend(('Pieper U', 'Webb B', 'Narayanan E', 'Sali A'))
        modpipe_software = ma.Software(
            name='ModPipe', classification='comparative modeling',
            location='https://salilab.org/modpipe/', type='program',
            version=self.modpipe_version,
            description='Comparative modeling pipeline')
        s.software.append(modpipe_software)
        modeller_software = ma.Software(
            name='MODELLER', classification='comparative modeling',
            location='https://salilab.org/modeller/', type='program',
            version=self.get_modeller_version(),
            citation=ihm.citations.modeller,
            description='Comparative modeling by satisfaction of '
                        'spatial restraints')
        s.software.append(modeller_software)

        tgtbeg = int(self.remarks['TARGET BEGIN'])
        tgtend = int(self.remarks['TARGET END'])
        if align:
            template_e = ma.Entity(align.template.primary,
                                   description="Template")
            s.entities.append(template_e)
            template_pdb = self.remarks['TEMPLATE PDB']
            # Some very old models use single-chain PDB templates with no
            # chain ID. These have all since been remediated, almost certainly
            # to use chain ID 'A'.
            template_chain = self.remarks['TEMPLATE CHAIN'] or 'A'
            tmpbeg, tmpend, tmpasym = self.get_mmcif_template_info(
                self.remarks['TEMPLATE BEGIN'], self.remarks['TEMPLATE END'],
                template_chain, template_pdb)
            template = ma.Template(entity=template_e, asym_id=tmpasym,
                                   model_num=1, name="Template structure",
                                   transformation=ma.Transformation.identity(),
                                   references=[ma.reference.PDB(template_pdb)])
        if align and align.target.primary == tgt_primary:
            target_e = template_e 
        else:
            target_e = ma.Entity(tgt_primary, description="Target")
            s.entities.append(target_e)
        target_e.references.extend(self.get_target_refs(tgtbeg, tgtend))
        chain_id = self.chain_id.strip() or 'A'
        asym = ma.AsymUnit(target_e, details='Model subunit', id=chain_id,
                           auth_seq_id_map=tgtbeg-1)
        s.asym_units.append(asym)
        asmb = ma.Assembly((asym,), name='Modeled assembly')

        if align:
            class OurAlignment(ma.alignment.Global, ma.alignment.Pairwise):
                pass
            p = ma.alignment.Pair(
                template=template.segment(align.template.gapped,
                                          tmpbeg, tmpend),
                target=asym.segment(align.target.gapped, tgtbeg, tgtend),
                score=ma.alignment.BLASTEValue(self.remarks['EVALUE']),
                identity=SequenceIdentity(self.remarks['SEQUENCE IDENTITY']))
            aln = OurAlignment(name="Modeling alignment",
                               software=modpipe_software, pairs=[p])
            s.alignments.append(aln)

        model = self.get_model_class(asym, self.atoms)(
                assembly=asmb, name='Best scoring model')
        model.qa_metrics.extend(self.get_scores(modeller_software,
                                                modpipe_software))

        model_group = ma.model.ModelGroup([model], name='All models')
        s.model_groups.append(model_group)

        protocol = ma.protocol.Protocol()
        mth = " %s" % align.template.method if align else ''
        protocol.steps.append(ma.protocol.TemplateSearchStep(
            name='ModPipe%s' % mth, software=modpipe_software,
            input_data=model, output_data=aln))
        protocol.steps.append(ma.protocol.ModelingStep(
            software=modeller_software, input_data=aln, output_data=model))
        protocol.steps.append(ma.protocol.ModelSelectionStep(
            software=modpipe_software, input_data=model, output_data=model))
        s.protocols.append(protocol)

        return s

    def get_target_refs(self, tgtbeg, tgtend):
        refmap = {'UniProt': ma.reference.UniProt,
                  'RefSeq': RefSeq,
                  'PlasmoDB': PlasmoDB}
        for db in self.seqdb:
            cls = refmap.get(db.name)
            if cls:
                yield cls(code=db.code, accession=db.accession,
                          align_begin=tgtbeg, align_end=tgtend)

    def get_model_class(self, asym, atoms):
        class MyModel(ma.model.Model):
            def get_atoms(self):
                pdb_resnum = None
                seqid = 1
                for a in atoms:
                    # Detect new residue if PDB resnum changed
                    pdb_this_resnum = a[22:26]
                    if pdb_resnum is not None and pdb_this_resnum != pdb_resnum:
                        seqid += 1
                    pdb_resnum = pdb_this_resnum
                    element = a[76:78].strip() or ihm.unknown
                    # Very old PDBs don't have sequence where the element should
                    # be, so ignore that; guess element as the first character
                    # of name
                    if a[73:76] == '1SG':
                        element = a[13:14].strip() or ihm.unknown
                    yield ma.model.Atom(asym_unit=asym, type_symbol=element,
                                        seq_id=seqid, atom_id=a[12:16].strip(),
                                        x=a[30:38].strip(), y=a[38:46].strip(),
                                        z=a[46:54].strip(),
                                        biso=a[60:66].strip(),
                                        occupancy=a[54:60].strip())
        return MyModel

    def get_scores(self, modeller_software, modpipe_software):
        tsvmod_method = self.remarks.get('TSVMOD METHOD')
        tsvmod_rmsd = self.remarks.get('TSVMOD RMSD')
        tsvmod_no35 = self.remarks.get('TSVMOD NO35')
        ga341 = self.remarks.get('GA341 SCORE',
                                 self.remarks.get('MODEL SCORE'))
        zdope = self.remarks.get('zDOPE SCORE',
                                 self.remarks.get('ZDOPE SCORE'))
        mpqs = self.remarks.get('MPQS',
                                self.remarks.get('MODPIPE QUALITY SCORE'))
        if not mpqs:
            return
        class MPQS(ma.qa_metric.Global, MPQSMetricType):
            name = "MPQS"
            description = "ModPipe Quality Score"
            software = modpipe_software

        class zDOPE(ma.qa_metric.Global, ma.qa_metric.ZScore):
            name = "zDOPE"
            description = "Normalized DOPE"
            software = modeller_software

        yield MPQS(mpqs)
        yield zDOPE(zdope)

        if tsvmod_rmsd:
            class TSVModRMSD(ma.qa_metric.Global, ma.qa_metric.Distance):
                name = "TSVMod RMSD"
                description = "TSVMod predicted RMSD (%s)" % tsvmod_method
                software = None

            class TSVModNO35(ma.qa_metric.Global, ma.qa_metric.NormalizedScore):
                name = "TSVMod NO35"
                description = ("TSVMod predicted native overlap (%s)"
                               % tsvmod_method)
                software = None

            yield TSVModRMSD(tsvmod_rmsd)
            yield TSVModNO35(tsvmod_no35)


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
    if args.pdb == '-':
        s = read_pdb(sys.stdin, r)
    else:
        with open(args.pdb) as fh:
            s = read_pdb(fh, r)
    system = s.get_system(args.align)

    if args.mmcif == '-':
        ma.dumper.write(sys.stdout, [system])
    else:
        if args.mmcif.endswith('.bcif'):
            mode, fmt = "wb", "BCIF"
        else:
            mode, fmt = "w", "mmCIF"
        with open(args.mmcif, mode) as fh:
            ma.dumper.write(fh, [system], format=fmt)
