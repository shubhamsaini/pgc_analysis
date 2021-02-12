#!/usr/bin/env python
"""

Author: Shubham Saini
shubhamsaini@eng.ucsd.edu

"""


import argparse
import sys
import pandas as pd

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--input", help="Input Summary Statistics File", required=True, type=str)
    parser.add_argument("--output", help="Filtered output", required=True, type=str)
    parser.add_argument("--avg-gp-file", help="File with AVG GP", required=False, type=str)
    parser.add_argument("--min-avg-gp", help="File with AVG GP", required=False, type=float, default=0.0)
    parser.add_argument("--min-af", help="File with AVG GP", required=False, type=float, default=0.0)
    args = parser.parse_args()

    MAF = args.min_af
    AVG_GP = args.min_avg_gp
    assoc_file = args.input
    output_file = args.output

    assoc_data = pd.read_csv(assoc_file, delim_whitespace=True, na_values="inf")
    assoc_data.dropna(inplace=True)

    if args.avg_gp_file:
        avg_gp_file = args.avg_gp_file
        avg_gp_data = pd.read_csv(avg_gp_file, delim_whitespace=True, names=['CHROM','BP','SNP','GP'])
        merged_data = pd.merge(assoc_data, avg_gp_data)
        merged_data = merged_data[merged_data['GP']>AVG_GP]
    else:
        merged_data = assoc_data[assoc_data['MAF']>MAF]

    merged_data.to_csv(output_file, sep="\t", columns=list(assoc_data.columns), index=False)

if __name__ == "__main__":
    main()
