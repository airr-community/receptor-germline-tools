# Utilities used in the other files in this package
from collections import namedtuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json


permitted_header_fields = [
    'release_version',  # Release version of the allele
    'allele_description_ref',  # Unique reference to the allele description, in standardized form (Repo:Label:Version)
    'aliases', # Alternative names for this sequence
    'locus', # Gene Locus (e.g. IGH)
    'chromosome', # Chromosome (e.g. 14)
    'sequence_type', # Sequence type (e.g. V)
    'functional', # Functional (T/F)
    'species', # Binomial designation of subject's species
    'species_subgroup', # Race, strain or other species subgroup to which this subject belongs
    'species_subgroup_type', # Type of species subgroup (e.g. strain)
    'subgroup_designation', # Identifier of the gene subgroup or clade, as (and if) defined
    'gene_designation', # Gene number or other identifier, as (and if) defined
    'allele_designation', # Allele number or other identifier, as (and if) defined
    'allele_similarity_cluster_designation', # ID of the similarity cluster used in this germline set, if designated
    'allele_similarity_cluster_member_id', # ID of the similarity cluster member used in this germline set, if designated
    'paralogs', # Labels of paralogs, if any
    'curational_tags', # Tags used to indicate curational status
    'gs_release_version', # Release version of the germline set
    'gs_germline_set_name', # Name of the germline set
    'gs_germline_set_ref', # Unique reference to the germline set, in standardized form (Repo:Label:Version)
]


AlleleData = namedtuple('AlleleData', 'name subgroup allele')

# is the name in IgLabel format?
def is_label(name):
    if name[4] != '-':       # ie subgroup specified
        return False

    label = name[5:].split('*')[0]

    if '-' in label or len(label) != 4:
        return False

    return True


# Return dummified name, given the call from the input file, or the current name, if it isn't in IgLabel format
def dummify(name, germline_data, dummy_subgroup, dummy_allele):
    if not is_label(name):
        return name

    if name in germline_data:
        if germline_data[name].subgroup:
            dummy_name = name[:4] + germline_data[name].subgroup + '-' + name[5:]
        else:
            dummy_name = name[:4] + dummy_subgroup + '-' + name[5:]

        if '*' not in name:
            dummy_name += '*' + dummy_allele
    else:
        dummy_name = name[:4] + dummy_subgroup + '-' + name[5:]

        if '*' not in name:
            dummy_name += '*' + dummy_allele

    return dummy_name


# Return undummified name, or the current name, if it isn't in dummified format
def undummify(name, dummy_allele):
    if '*' not in name:
        return name

    if '-' not in name:
        return name

    gene_number = name.split('*')[0]
    gene_number = gene_number[name.index('-') + 1:]

    if len(gene_number) != 4 or '-' in gene_number:
        return name

    label_form = name[:4] + '-' + gene_number

    allele = name.split('*')[1:]
    allele = ''.join(allele)

    if allele == dummy_allele:
        undummy_name = label_form
    else:
        undummy_name = label_form + '*' + allele

    return undummy_name


# read fasta into dict
def read_fasta(infile):
    res = {}
    recs = SeqIO.parse(infile, 'fasta')
    for rec in recs:
        res[rec.id] = str(rec.seq).upper()
    return res


# convert our dict format to SeqRecords
def toSeqRecords(seqs):
    recs = []
    for name, seq in seqs.items():
        recs.append(SeqRecord(Seq(seq), id=name, description=''))
    return recs


# write fasta from dict
def write_fasta(outfile, seqs):
    seqs = {k: v for k, v in sorted(seqs.items(), key=lambda item: item[0])}
    SeqIO.write(toSeqRecords(seqs), outfile, 'fasta')


# Read germline data we want into AlleleDescriptions
def read_germline_data(args):
    germline_data = {}

    with open(args.germline_set, 'r') as fi:
        germline_set = json.load(fi)

    for allele_description in germline_set['GermlineSet']['allele_descriptions']:
        subgroup = allele_description['subgroup_designation'] if (allele_description['subgroup_designation']) else args.dummy_subgroup
        allele = allele_description['allele_designation'] if allele_description['allele_designation'] else args.dummy_allele
        germline_data[allele_description['label']] = AlleleData(allele_description['label'], subgroup, allele)

    return germline_data

