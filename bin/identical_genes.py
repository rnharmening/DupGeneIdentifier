#!/usr/bin/env python3

import pandas as pd
import argparse
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Find identical genes from a duplicates.csv',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

    parser.add_argument("-i", "--input", nargs="+", type=str,
        help="A list of *duplicates.tsv file, seperated by an space")
    parser.add_argument("-o", "--output", nargs="+", type=str,
        help="List of output files, has to be the same lenght as the " \
             "input file list. The index in each list will be matched")
    
    args = parser.parse_args()
    
    if args.output is None or args.input is None:
        parser.error("-i and -o have to be specified (smae number of elements)")

    if len(args.input) != len(args.output):
        parser.error('The number of elements passed to --input and --output has to be the same.')
    
    return args


if __name__ == "__main__":
    args = parse_arguments()
    in_files = args.input
    out_files = args.output

    to_delete, df = None, None

    for in_, out_ in zip(in_files, out_files):

        # Load and pre filter dataframe
        df = pd.read_csv(in_, sep="\t")
        df = df[df["pident"] == 100]
        df = df[df["qcovs"] == 100]
        df = df[df["qaccver"] != df["saccver"]]

        df = df.sort_values(by=["qaccver", "saccver"])

        # get the set of genes for removal, due to them being allready coverred by an identical gene
        q_s_set = set(zip(df.qaccver, df.saccver))
        s_q_set = set(zip(df.saccver, df.qaccver))
        to_delete = pd.DataFrame({"toDelete": list(set(max(tpl) for tpl in q_s_set.intersection(s_q_set)))})

        to_delete.to_csv(out_, sep="\n", index=False)

        