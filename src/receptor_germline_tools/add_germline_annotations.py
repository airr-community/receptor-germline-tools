#!python
# Add germline annotations (subgroup, allele etc) to an alignment file, using metadata from an AIRRC germline set

# Copyright (c) 2022 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE


import csv
import argparse


def main():
    parser = argparse.ArgumentParser(description='Add germline annotations to an alignment file, using metadata from an AIRRC germline set')
    parser.add_argument('input_file', help='alignments to annotate (csv, tsv)')
    parser.add_argument('output_file', help='output with added annotations')
    parser.add_argument('germline_set', help='AIRR standard germline set to use for metadata (JSON)')
    args = parser.parse_args()
    
