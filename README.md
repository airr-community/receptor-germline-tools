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

Two command-line tools are provided. `add_germline_annotations` is intended to operate on the output of a sequence annotation tool,
but is sufficiently flexible to operate on most csv or tsv files that have columns containing
allele names. It will convert such names in nominated columns between label format and the
'dummy' format containing subgroup and allele. `convert_fasta_labels` will perform the
same operations on sequence names in a FASTA file. Further details are provided below.

By default, the dummy values used by the tools are 0 for the subgroup, and 00 for the allele. These values can be changed if
necessary. If you do need to change the values, please use values that are not used by any existing allele in the
reference set, so that the dummy values can be distinguished.

# Installation
Python v3.9 or greater is required.

```bash
pip install biopython   # if not already installed
pip install receptor-germline-tools
```


### add_germline_annotations

```commandline
usage: add_germline_annotations [-h] [-c CALL_COLUMNS] [-s DUMMY_SUBGROUP] [-a DUMMY_ALLELE] [-u] input_file output_file germline_set

Convert IgLabel-style labels in an alignment file to dummy IUIS format

positional arguments:
  input_file            alignments to annotate (csv, tsv)
  output_file           output with added annotations
  germline_set          AIRR standard germline set to use for metadata (JSON)

optional arguments:
  -h, --help            show this help message and exit
  -c CALL_COLUMNS, --call_columns CALL_COLUMNS
                        Names of one or more columns to be processed, separated by commas (default: v_call,d_call,j_call)
  -s DUMMY_SUBGROUP, --dummy_subgroup DUMMY_SUBGROUP
                        The subgroup to be used, where no value is specified in the germline set metadata (default: 0)
  -a DUMMY_ALLELE, --dummy_allele DUMMY_ALLELE
                        The allele to be used, where no value is specified in the germline set metadata (default: 00)
  -u, --un_dummy        translate dummy names back to label form
```

#### Description:

The input file should be comma or tab-separated data, in which one or more columns contain allele
names. The format of the file is determined automatically from its content.
MiAIRR column names are used by default. In the absence of the -u option, any allele names that match temporary labels that
are listed in the germline set will be converted into a 'dummy IUIS format',
including subgroup and allele number. These numbers are taken from the germline set metadata. Names that do not match the temporary
label format will be left as-is. Translated names (and those left unchanged) are put in a new column with the suffix '_dummy'.
If the metadata does not provide values, the 'dummy' values are used: by default these are 0 for subgroup, and 00 for allele. If the -u option is
specified, the operation is reversed: any names that match the dummy label format are
converted back to the 'pure' label format, without subgroup or allele.

Cells in the columns to be processed may contain multiple allele names. In this case the names 
should be separated by a comma.

### convert_fasta_labels

```commandline
usage: convert_fasta_labels [-h] [-s DUMMY_SUBGROUP] [-a DUMMY_ALLELE] [-u] input_file output_file germline_set

Convert IgLabel-style labels in a FASTA file to dummy IUIS format

positional arguments:
  input_file            records to convert (FASTA)
  output_file           converted output (FASTA)

optional arguments:
  -h, --help            show this help message and exit
  -g GERMLINE_SET, --germline_set GERMLINE_SET 
                        AIRR standard germline set to use for metadata (JSON)
  -s DUMMY_SUBGROUP, --dummy_subgroup DUMMY_SUBGROUP
                        subgroup to use when no subgroup has been defined (default: 0)
  -a DUMMY_ALLELE, --dummy_allele DUMMY_ALLELE
                        allele to use when no allele has been defined (default: 00)
  -u, --un_dummy        translate dummy names back to label form
```

#### Description:

The input file should be a FASTA file. In the absence of the -u option, any sequence names that match temporary labels that
are listed in the germline set will be converted into a 'dummy IUIS format', including subgroup and allele number. These numbers are taken from the germline set metadata.
If the metadata does not provide values, the 'dummy' values are used: by default these are 0 for subgroup, and 00 for allele. If the -u option is
specified, the operation is reversed: any names that match the dummy label format are
converted back to the 'pure' label format, without subgroup or allele.
