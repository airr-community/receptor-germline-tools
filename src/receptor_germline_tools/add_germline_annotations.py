#!python
# Add germline annotations (subgroup, allele etc) to an alignment file, using metadata from an AIRRC germline set

# Copyright (c) 2022 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE


import csv
import argparse
import os
import json
from collections import namedtuple
import csv

AlleleData = namedtuple('AlleleData', 'name subgroup allele')

# is the lanem in IgLabel format?
def is_label(name):
    if name[4] != '-':       # ie subgroup specified
        return False

    label = name[5:].split('*')[0]

    if '-' in label or len(label) != 4:
        return False

    return True


# Return dummified name, given the call from the input file
def dummify(name, germline_data, dummy_subgroup, dummy_allele):
    if is_label(name):
        if name in germline_data:
            dummy_name = name[:4] + germline_data[name].subgroup + '-' + name[5:]

            if '*' not in name:
                dummy_name += '*' + dummy_allele

            return dummy_name
        else:
            dummy_name = name[:5] + dummy_subgroup + name[5:]
            return dummy_name
    else:
        return name


def main():
    parser = argparse.ArgumentParser(description='Add germline annotations to an alignment file, using metadata from an AIRRC germline set')
    parser.add_argument('input_file', help='alignments to annotate (csv, tsv)')
    parser.add_argument('output_file', help='output with added annotations')
    parser.add_argument('germline_set', help='AIRR standard germline set to use for metadata (JSON)')
    parser.add_argument('-v', '--v_call', default='v_call')
    parser.add_argument('-d', '--d_call', default='d_call')
    parser.add_argument('-j', '--j_call', default='j_call')
    parser.add_argument('-s', '--dummy_subgroup', default='0')
    parser.add_argument('-a', '--dummy_allele', default='00')

    args = parser.parse_args()

    missing_files = False
    for filespec in [(args.input_file, 'input_file'), (args.germline_set, 'germline_set')]:
        if not os.path.isfile(filespec[0]):
            print(f'{filespec[1]} {filespec[0]} does not exist.')
            missing_files = True

    if missing_files:
        exit(1)

    germline_data = {}

    with open(args.germline_set, 'r') as fi:
        germline_set = json.load(fi)

    for allele_description in germline_set['GermlineSet']['allele_descriptions']:
        # A blank subgroup designation is used for example in J genes were no subgroup is ever assigned
        # A null subgroup indicates that a subgroup designation has not been determined
        subgroup = allele_description['subgroup_designation'] if (allele_description['subgroup_designation'] or allele_description['subgroup_designation'] == "") else args.dummy_subgroup
        allele = allele_description['allele_designation'] if allele_description['allele_designation'] else args.dummy_allele
        germline_data[allele_description['label']] = AlleleData(allele_description['label'], subgroup, allele)

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
                for call in (args.v_call, args.d_call, args.j_call):
                    if call in headers:
                        headers.append(call + '_dummy')

                writer = csv.DictWriter(fo, fieldnames=headers)
                writer.writeheader()

            for call in (args.v_call, args.d_call, args.j_call):
                if call in row:
                    names = [dummify(name.replace(' ', ''), germline_data, args.dummy_subgroup, args.dummy_allele) for name in row[call].split(',')]
                    row[call + '_dummy'] = ','.join([name for name in names])
                    print(f"{row[call]} -> {row[call + '_dummy']}")


if __name__ == "__main__":
    main()
