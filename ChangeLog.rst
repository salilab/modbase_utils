1.0 - 2024-09-30
================
 - Information about the database sequence for the model (e.g. in
   UniProt) is now provided both in the ModelCIF-specific
   ``ma_target_ref_db_details`` table and in the more general
   ``struct_ref`` table.

0.5 - 2024-01-12
================
 - New ModBase models containing multiline TITLE records and new-style
   EXPDTA records (where Modeller version information has been moved
   to a REMARK) are now handled correctly.

0.4 - 2022-11-16
================
 - Author names are now PDB style ("Lastname, A.B.") rather than
   PubMed style ("Lastname AB").
 - Old ModBase models containing UNK residues are now handled correctly.

0.3 - 2022-05-10
================
 - The ModPipe databases (PDB95, UniProt90, profiles) are now linked from
   the modeling protocol in the output mmCIF file.

0.2 - 2022-04-15
================
 - Sequence information for templates is now only written to template-specific
   categories in the mmCIF/BinaryCIF, not to the entity, entity_poly etc.
   tables, to properly comply with the ModelCIF dictionary.

0.1 - 2022-03-29
================
 - First stable release. ModBase PDB files, together with the modeling
   alignments, can be converted to mmCIF compliant with the latest version
   (1.3.6) of the ModelCIF dictionary.
