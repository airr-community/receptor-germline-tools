.. _working_with_labels:

Working with IgLAbel-style labels
=================================

Germline sets published on `OGRDB <https://ogrdb.airr-community.org/>`_ may contain sequences which have not been allocated official names 
by `IUIS <https://iuis.org/committees/nom/immunoglobulins-ig-t-cell-receptors-tr-and-major-histocompatibility-mh-nomenclature-sub-committee/>`_. 
These sequences will normally be allocated a four letter label, as outlined in `Lees et al., 2023 <https://www.sciencedirect.com/science/article/pii/S2667119023000058>`_.
In some cases, these sequences may not have been mapped to a specific gene. This can happen, for example, where the sequence
has been discovered by inference from AIRR-seq repertoires. The mouse sets on OGRDB are a good example of this, as they are derived
largely from AIRR-seq data. The human sets follow a more traditional naming scheme in order to remain as aligned as possible with the IUIS nomenclature.

Because some existing tools have difficulty working with allele names that do not include a subgroup number and/or allele number, in FASTA sets provided by OGRDB,
a dummy subgroup and/or allele number is provided so that the sequences can be processed. In the mouse sets, for example, names of the form IGKV0-2GER\*00 are used,
where the 0 and 00 represent the dummy values. In the AIRR-C JSON format, the label would be shown as IGKV-2GER, and the subgroup designation and allele designation
would be null. Because the dummy values are not part of the name, and are provided only to facilitate toolchains, we strongly recommend that they are removed before publication,
so that the published names are always consistent and follow the temporary label format.

In case they may be useful, we provide two tools here which can convert between the dummy format and the temporary label format.

:ref:`add_germline_annotations` is intended to operate on the output of a sequence annotation tool, but is sufficiently flexible to operate on most csv or tsv files that have columns containing
allele names. It will convert such names in nominated columns between label format and the 'dummy' format containing subgroup and allele. 

:ref:`convert_fasta_labels` will perform the same operations on sequence names in a FASTA file. 

Both tools will operate on a file containing a mixture of dummy and real values, and will only convert the dummies, provided care is taken to pick
values for dummies that are not otherwise used. The defaults of 0 and 00 are generally safe to use, as they are never used in IUIS-allocated names.
