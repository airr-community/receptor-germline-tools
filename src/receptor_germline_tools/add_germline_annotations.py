#!python
# Convert IgLabel-style labels to 'dummy' IUIS format

# Copyright (c) 2022 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE


import csv
import argparse
import os
import csv
from receptor_germline_tools.germline_utils import *


def main():
    parser = argparse.ArgumentParser(description='Convert IgLabel-style labels in an alignment file to dummy IUIS format')
    parser.add_argument('input_file', help='alignments to annotate (csv, tsv)')
    parser.add_argument('output_file', help='output with added annotations')
    parser.add_argument('germline_set', help='AIRR standard germline set to use for metadata (JSON)')
    parser.add_argument('-c', '--call_columns', help='Names of one or more columns to be processed, separated by commas', default='v_call,d_call,j_call')
    parser.add_argument('-s', '--dummy_subgroup', help='The subgroup to be used, where no value is specified in the germline set metadata', default='0')
    parser.add_argument('-a', '--dummy_allele', help='The allele to be used, where no value is specified in the germline set metadata', default='00')
    parser.add_argument('-u', '--un_dummy', help='translate dummy names back to label form', action="store_true")

    args = parser.parse_args()
    call_columns = args.call_columns.split(',')

    missing_files = False
    for filespec in [(args.input_file, 'input_file'), (args.germline_set, 'germline_set')]:
        if not os.path.isfile(filespec[0]):
            print(f'{filespec[1]} {filespec[0]} does not exist.')
            missing_files = True

    if missing_files:
        exit(1)

    germline_data = read_germline_data(args)

    with open(args.input_file, 'r') as fi, open(args.output_file, 'w', newline='') as fo:
        writer = None
        tabfile = '\t' in fi.readline()
        fi.seek(0)

        if tabfile:
            reader = csv.DictReader(fi, delimiter='\t')
        else:
            reader = csv.DictReader(fi)

        for row in reader:
            if not writer:
                headers = list(reader.fieldnames)
                for call in call_columns:
                    if call in headers:
                        if args.un_dummy:
                            headers.append(call + '_undummy')
                        else:
                            headers.append(call + '_dummy')

                writer = csv.DictWriter(fo, fieldnames=headers)
                writer.writeheader()

            for call in call_columns:
                if call in row:
                    if args.un_dummy:
                        names = [undummify(name.replace(' ', ''), args.dummy_allele) for name in row[call].split(',')]
                    else:
                        names = [dummify(name.replace(' ', ''), germline_data, args.dummy_subgroup, args.dummy_allele) for name in row[call].split(',')]

                    if args.un_dummy:
                        row[call + '_undummy'] = ','.join([name for name in names])
                        # print(f"{row[call]} -> {row[call + '_undummy']}")
                    else:
                        row[call + '_dummy'] = ','.join([name for name in names])
                        # print(f"{row[call]} -> {row[call + '_dummy']}")

            writer.writerow(row)



if __name__ == "__main__":
    main()
