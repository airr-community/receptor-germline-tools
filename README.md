# Tools for working with AIRR standard germline sets

AIRR standard germline sets (such as those published on [OGRDB](https://ogrdb.airr-community.org/)) contain rich metadata for each included allele. This metadata may not always be provided in the 
name. In summary:

- Germline sets may contain sequences with names that have been allocated by IUIS. These are 'IMGT style' names and will match those used by IMGT.
- Germline sets may also contain sequences that have not yet been allocated names by IUIS. These are provided with 'temporary labels' in a format that
is designed to be consistent, and readily distinguishable from names allocated by IUIS. They are based around a four-character identifier, for example IGHV-H7DF.
- Where two or more alleles of the same gene have been identified, the temporary label will include an allele number, for example IGHV-H7DF&ast;01, IGHV-H7DF&ast;02.
- A subgroup (family) number is never included in the temporary label.

For further details on the naming convention used, and the rationale behind it, please refer to our [poster](https://wordpress.vdjbase.org/index.php/ogrdb_news/germline-set-creation-and-naming/).

We encourage authors of tools and those building analysis pipelines to use the metadata provided in the germline set, rather than parsing the allele name to find the allele number, gene number and so on. This is simpler
programatically, and also allows for the possibility that some information is not available, for reasons outlined in the poster referenced above. Likewise we encourage
authors of tools to avoid assumptions about the format of a name, and in particular to accept the temporary label format.

Because some existing tools may have difficulty working with allele names that do not include a subgroup number and/or allele number, we provide here tools that will convert allele names into a format
that should be acceptable to such tools by inserting dummy values. We hope that these will be less needed over time. The tools can also remove dummy values. We strongly encourage their removal before 
publication, so that the published names are always consistent and follow the temporary label format.

Two tools are provided. `add_germline_annotations` is intended to operate on the output of a sequence annotation tool,
but is sufficiently flexible to operate on most csv or tsv files that have columns containing
allele names. It will convert such names in nominated columns between label format and the
'dummy' format containing subgroup and allele. `convert_fasta_labels` will perform the
same operations on sequence names in a FASTA file. Further details are provided below.

# Installation

```bash
pip install receptor-germline-tools
```
The module requires [Biopython](https://biopython.org).


### add_germline_annotations

#### Usage:

#### Description:

The input file should be comma or tab-separated data, in which one or more columns contain allele
names. MiAIRR tsv has been tested and its `v/d/j_call` column names are used by default. In the absence of the -u option, any allele names that match temporary labels that
are listed in the germline set will be converted into a 'dummy IUIS format',
including subgroup and allele number. These numbers are taken from the germline set metadata.
If the metadata does not provide values, the 'dummy' values are used. If the -u option is
specified, the operation is reversed: any names that match the dummy label format are
converted back to the 'pure' label format, without subgroup or allele.

Cells in the columns to be processed may contain multiple allele names. In this case the names 
should be separated by a comma.

### convert_fasta_labels

#### Usage:

#### Description:

The input file should be a FASTA file. In the absence of the -u option, any sequence names that match temporary labels that
are listed in the germline set will be converted into a 'dummy IUIS format', including subgroup and allele number. These numbers are taken from the germline set metadata.
If the metadata does not provide values, the 'dummy' values are used. If the -u option is
specified, the operation is reversed: any names that match the dummy label format are
converted back to the 'pure' label format, without subgroup or allele.
