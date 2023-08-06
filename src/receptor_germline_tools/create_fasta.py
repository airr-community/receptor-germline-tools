#!python
# Create a FASTA file from an AIRR-C formatted germline set

# Copyright (c) 2023 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE

import argparse
import json
from receptor_germline_tools.germline_utils import write_fasta, permitted_header_fields



def ornull(rec, field):
    if field in rec and rec[field] is not None:
        return rec[field]
    else:
        return ''
    
    
def append_if_present(result, rec, field, force, prefix=''):
    val = ornull(rec, field)
    if val != '' or force:
        try:
            result.append(f"{prefix+field}={','.join(val)}")
        except TypeError:
            result.append(f"{prefix+field}={val}")
    

def main():
    parser = argparse.ArgumentParser(description='Convert IgLabel-style labels in a FASTA file to dummy IUIS format')
    parser.add_argument('germline_set', help='Germline set in AIRR-C format')
    parser.add_argument('output_file', help='Germline set in FASTA format')
    parser.add_argument('-g', '--gapped_set', help='Create IMGT-gapped set', action="store_true", default=False)
    parser.add_argument('-f', '--fields', help='add comma-separated fields to headers, e.g. release_version,aliases')
    parser.add_argument('-n', '--include_null', help='include fields with null values', action="store_true", default=False)

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

    gs_params = []

    for param in permitted_header_fields:
        if param.startswith('gs_') and param in required_header_fields:
            append_if_present(gs_params, germline_set, param[3:], args.include_null, prefix='gs_')

    seqs = {}

    for allele_description in germline_set['allele_descriptions']:
        header = [allele_description['label'] if 'label' in allele_description else allele_description['id'] if 'id' in allele_description else 'XXX']
        header.extend(gs_params)

        for param in required_header_fields:
            if not param.startswith('gs_'):
                append_if_present(header, allele_description, param, args.include_null)
        
        seq = None
        if args.gapped_set and 'sequence_type' in allele_description \
            and allele_description['sequence_type'] == 'V' \
            and 'v_gene_delineations' in allele_description:
            for delineation in allele_description['v_gene_delineations']:
                if 'delineation_scheme' in delineation and delineation['delineation_scheme'] == 'IMGT' and 'aligned_sequence' in delineation:
                    seq = delineation['aligned_sequence']
                    break

            if seq is None:
                print(f"No IMGT-gapped sequence is provided in the germline set for {header[0]}")
                continue

        if seq is None:
            if 'coding_sequence' in allele_description:
                seq = allele_description['coding_sequence']
            else:
                print(f"No coding sequence is provided in the germline set for {header[0]}")
                continue

        seqs[' | '.join(header)] = seq

    
    write_fasta(args.output_file, seqs)


if __name__ == "__main__":
    main()
