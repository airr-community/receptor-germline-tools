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
programatically, and also allows for the possibility that some information is not available, for reasons outlined in the poster referenced above.
To support this, we provide here a command-line tool, add_germline_annotations, which will add commonly-used metadata to a rearrangement file. The tool is intended to be used with AIRR standard rearrangement files, 
but is sufficientky flexible to work with most rearrangement schemas in CSV or TSV format.

Because some existing tools may have difficulty working with allele names that do not include a subgroup number and/or allele number, we also provide tools that will convert allele names into a format
that should be acceptable to such tools by inserting dummy values. We hope that these will be less needed over time. The tools can also remove dummy values. We strongly encourage their removal before 
publication, so that the published names are always consistent and follow the temporary label format.

# Installation

```bash
pip install receptor-germline-tools
```

The module requires [Biopython](https://biopython.org).


### add_germline_annotations

