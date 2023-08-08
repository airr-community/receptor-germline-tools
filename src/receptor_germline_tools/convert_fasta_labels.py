#!python
# Convert IgLabel-style labels in a FASTA file to dummy IUIS format

# Copyright (c) 2022 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE


import csv
import argparse
import os
import json
from collections import namedtuple
import csv

try:
    from receptor_germline_tools.germline_utils import *
except:
    from germline_utils import *


def get_parser():
    parser = argparse.ArgumentParser(description='Convert IgLabel-style labels in a FASTA file to dummy IUIS format')
    parser.add_argument('input_file', help='records to convert (FASTA)')
    parser.add_argument('output_file', help='converted output (FASTA)')
    parser.add_argument('-g', '--germline_set', help='AIRR standard germline set to use for metadata (JSON)', default=None)
    parser.add_argument('-s', '--dummy_subgroup', help='subgroup to use when no subgroup has been defined', default='0')
    parser.add_argument('-a', '--dummy_allele', help='allele to use when no allele has been defined', default='00')
    parser.add_argument('-u', '--un_dummy', help='translate dummy names back to label form', action="store_true", default=False)
    return parser


def main():
    args = get_parser().parse_args()

    if not os.path.isfile(args.input_file):
        print(f'{args.input_file} does not exist.')
        exit(0)

    germline_data = {}
    if args.germline_set:
        if os.path.isfile(args.germline_set):
            germline_data = read_germline_data(args)
        else:
            print(f'{args.germline_set} does not exist.')
            exit(0)

    recs = read_fasta(args.input_file)
    converted_recs = {}

    for name, seq in recs.items():
        if '|' in name:
            name = name.split('|')[0].replace(' ', '')
        if args.un_dummy:
            converted_recs[undummify(name, args.dummy_allele)] = seq
        else:
            converted_recs[dummify(name, germline_data, args.dummy_subgroup, args.dummy_allele)] = seq

    write_fasta(args.output_file, converted_recs)


if __name__ == "__main__":
    main()
