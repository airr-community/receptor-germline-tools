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
from germline_utils import *


def main():
    parser = argparse.ArgumentParser(description='Convert IgLabel-style labels in a FASTA file to dummy IUIS format')
    parser.add_argument('input_file', help='records to convert (FASTA)')
    parser.add_argument('output_file', help='converted output (FASTA)')
    parser.add_argument('germline_set', help='AIRR standard germline set to use for metadata (JSON)')
    parser.add_argument('-s', '--dummy_subgroup', default='0')
    parser.add_argument('-a', '--dummy_allele', default='00')
    parser.add_argument('-u', '--un_dummy', help='translate dummy names back to label form', action="store_true")

    args = parser.parse_args()

    missing_files = False
    for filespec in [(args.input_file, 'input_file'), (args.germline_set, 'germline_set')]:
        if not os.path.isfile(filespec[0]):
            print(f'{filespec[1]} {filespec[0]} does not exist.')
            missing_files = True

    if missing_files:
        exit(1)

    germline_data = read_germline_data(args)

    recs = read_fasta(args.input_file)
    converted_recs = {}

    for name, seq in recs.items():
        if '|' in name:
            name = name.split('|')[0].replace(' ', '')
        if args.un_dummy:
            converted_recs[undummify(name, args.dummy_allele)] = seq
        else:
            converted_recs[dummify(name, germline_data, args.dummy_subgroup, args.dummy_allele)] = seq

    write_fasta(converted_recs, args.output_file)


if __name__ == "__main__":
    main()
