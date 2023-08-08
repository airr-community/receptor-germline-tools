Introduction
============

This package contains utilities which may be useful in working with AIRR-C JSON germline sets

Installation
------------

The package requires python, and can be installed using pip:

``pip install receptor-germline-tools``

It requires `Biopython <https://biopython.org/>`_, which may need to be installed separately.

Using the rich information in AIRR-C germline sets
--------------------------------------------------

Quick summary:

Given an AIRR-C germline set ``germlines.json``, create a FASTA-formatted IMGT-gapped set, in which the header lists paralogs (other items
in the set with the same sequence) and aliases (other names by which this item has been known)::

   create_fasta -i germlines.json -o germlines.fasta -g -f paralogs,aliases

Given an AIRR-C rearrangements file or Changeo file ``rearrangements.tsv``, create a copy that also includes columns for subgroup number, 
gene number and allele number of each V,D and J allele, using information from the germline set::

   add_germline_annotations rearrangements.tsv germlines.json \ 
     rearrangements_annotated.tsv \ 
     -f subgroup_designation,gene_designation,allele_designation


:ref:`list_of_fields`

More detail:

:ref:`create_fasta` - Create a gapped or ungapped FASTA file from an AIRR-C JSON-format germline set, optionally containing additional information in the header

:ref:`annotate_rearrangements` - Annotate an AIRR-C rearrangements file with additional germline information from a AIRR-C JSON-format germline set

Working with IGLabel-style sequence names
-----------------------------------------

Alleles in the mouse germline sets are allocated four-letter  `IgLabel-style <https://github.com/williamdlees/IgLabel>`_ labels with 'dummy' subgroup and allele designations.
These can be used as-is with many tools.

:ref:`Further information<working_with_labels>`

:ref:`add_germline_annotations` - Convert IgLabel-style labels in an alignment file to dummy IUIS format

:ref:`convert_fasta_labels` - Create a custom auxiliary data file for IgBlast

