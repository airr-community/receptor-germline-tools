# Utilities used in the other files in this package
from collections import namedtuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json


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
            dummy_name = name[:4] + germline_data[name].subgroup + name[5:]

        if '*' not in name:
            dummy_name += '*' + dummy_allele

    else:
        dummy_name = name

    return dummy_name


# Return undummified name, or the current name, if it isn't in dummified format
def undummify(name, dummy_allele):
    if '*' not in name:
        return name

    if name[3] == 'J' and '-' not in name:
        gene_number = name.split('*')[0]
        gene_number = gene_number[4:]
    else:
        if '-' not in name:
            return name

        gene_number = name.split('*')[0]
        gene_number = gene_number[name.index('-') + 1:]

    if len(gene_number) != 4:
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
def write_fasta(seqs, outfile):
    SeqIO.write(toSeqRecords(seqs), outfile, 'fasta')


# Read germline data we want into AlleleDescriptions
def read_germline_data(args):
    germline_data = {}

    with open(args.germline_set, 'r') as fi:
        germline_set = json.load(fi)

    for allele_description in germline_set['GermlineSet']['allele_descriptions']:
        if allele_description['sequence_type'] != 'J':
            subgroup = allele_description['subgroup_designation'] if (allele_description['subgroup_designation']) else args.dummy_subgroup
        else:
            subgroup = ''
        allele = allele_description['allele_designation'] if allele_description['allele_designation'] else args.dummy_allele
        germline_data[allele_description['label']] = AlleleData(allele_description['label'], subgroup, allele)

    return germline_data

