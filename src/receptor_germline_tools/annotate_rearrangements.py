#!python
# Annotate an AIRR-standard rearrangements file with selected metadata from an AIRR-standard germline set

# Copyright (c) 2023 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE

import argparse
import json
import csv

try:
    from receptor_germline_tools.germline_utils import permitted_header_fields
except:
    from germline_utils import permitted_header_fields


def main():
    parser = argparse.ArgumentParser(description='Convert IgLabel-style labels in a FASTA file to dummy IUIS format')
    parser.add_argument('input_file', help='Annotation file in TSV format')
    parser.add_argument('germline_set', help='Germline set in AIRR-C format')
    parser.add_argument('output_file', help='Extended annotation file in TSV format')
    parser.add_argument('-f', '--fields', help='add comma-separated fields to headers, e.g. release_version,aliases')
    parser.add_argument('-d', '--distinct', help='where there are multiple calls, provide distinct values for each call, semicolon separated', action="store_true", default=False)

    args = parser.parse_args()
    required_header_fields = []
    if args.fields:
        required_header_fields = args.fields.split(',')

    for field in required_header_fields:
        if field not in permitted_header_fields:
            print(f'Unrecognized header field: {field}')
            exit(1)

    with open(args.germline_set, 'r') as fi:
        germline_set = json.load(fi)

    if len(germline_set['GermlineSet']) > 0:
        germline_set = germline_set['GermlineSet'][0]
    else:
        germline_set = germline_set['GermlineSet']

    gs_params = {}

    for param in permitted_header_fields:
        if param.startswith('gs_') and param in required_header_fields:
            gs_params[param] = germline_set[param[3:]]

    annot_lookup = make_annot_lookup(germline_set, required_header_fields)

    with open(args.input_file, 'r', newline='\n') as fi, open(args.output_file, 'w', newline='\n') as fo:
        reader = csv.DictReader(fi, delimiter='\t')
        fieldnames = reader.fieldnames
        fieldnames.extend(gs_params)

        for seq_type in 'v', 'd', 'j':
            if seq_type + '_call' in fieldnames:
                for field in required_header_fields:
                    if not field.startswith('gs_'):
                        fieldnames.append(f"{seq_type}_{field}")

        writer = csv.DictWriter(fo, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in reader:
            row.update(gs_params)

            for seq_type in 'v', 'd', 'j':
                if seq_type + '_call' in row:
                    for field in required_header_fields:
                        if not field.startswith('gs_'):
                            val = []
                            for call in row[seq_type + '_call'].replace(' ', '').split(','):
                                if call in annot_lookup and field in annot_lookup[call]:
                                    val.append(annot_lookup[call][field])
                            if args.distinct:
                                row[f"{seq_type}_{field}"] = ';'.join(val)
                            else:
                                val = list(set(val))
                                row[f"{seq_type}_{field}"] = ','.join([ x for x in val if x])

                    
            writer.writerow(row)


def make_annot_lookup(germline_set, required_header_fields):
    annot_lookup = {}

    for allele_description in germline_set['allele_descriptions']:
        label = allele_description['label'] if 'label' in allele_description else allele_description['id'] if 'id' in allele_description else None

        if not label:
            continue
        
        annot_lookup[label] = {}

        for field in required_header_fields:
            if not field.startswith('gs_'):
                if field in allele_description:
                    try:
                        annot_lookup[label][field] = ','.join([x for x in allele_description[field] if x])
                       
                    except:
                        annot_lookup[label][field] = str(allele_description[field]) if allele_description[field] else ''

    return annot_lookup


if __name__ == "__main__":
    main()

