#!/usr/bin/env python3

##
#
# Author: Robin Harmening
#
##

import pandas as pd
import argparse
from pathlib import Path
import os


def filter_blast_df(df: pd.DataFrame, *, pident = 80, qcovs = 75):

    df.columns = "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen qcovs qcovhsp".split()

    #remove HSPs where query and subject id are identical
    filtered = df[df['qaccver'] != df['saccver']]

    filtered = filtered[filtered.pident > pident]
    filtered = filtered[filtered.qcovs > qcovs]

    return filtered


def get_base_name_no_ext(path: str):
    return os.path.splitext(os.path.basename(path))[0]


def parse_args():
    parser = argparse.ArgumentParser(
    description='Find Duplicated Genes based on blast output table.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("-i", "--input", nargs="+", type=str,
        help="A list of blast output table files (.tsv), where the columns are in this order:"\
            "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen qcovs qcovhsp")
    parser.add_argument("-o", "--output", nargs="+", type=str, 
        help="List of output file names, has to be same length as input list")
    parser.add_argument("--pident", type=int, default=80,
        help="The threshold 'pident' value to consider two genes as duplicates (>pident)")
    parser.add_argument("--qcovs", type=int, default=75,
        help="The threshold 'qcovs' value to consider two genes as duplicates (>qcovs)")

    args = parser.parse_args()
    
    if args.output is None or args.input is None:
        parser.error("-i and -o have to be specified (smae number of elements)")

    if len(args.input) != len(args.output):
        parser.error('The number of elements passed to --input and --output has to be the same.')
    
    return args


if __name__ == "__main__":
    
    args = parse_args()

    for in_file, out_file in zip(args.input, args.output):
        df = pd.read_csv(in_file, header=None, sep="\t")
        filtered = filter_blast_df(df, pident=args.pident, qcovs=args.qcovs)
        filtered.to_csv(out_file, index=False, sep="\t")
