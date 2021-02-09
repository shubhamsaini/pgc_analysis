#!/usr/bin/env python
"""

Author: Shubham Saini
shubhamsaini@eng.ucsd.edu

"""


import argparse
import sys
import numpy as np
from cyvcf2 import VCF

def gp_position(j, k):
    p1 = (k*(k+1)/2)+j
    p2 = (j*(j+1)/2)+k
    return int(max(p1, p2))

def gp_position_np(row):
    j,k = row[0], row[1]
    p1 = (k*(k+1)/2)+j
    p2 = (j*(j+1)/2)+k
    return int(max(p1, p2))

def PrintLine(text, f):
    f.write(text+"\n")
    f.flush()

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--vcf", help="Argument 1 description", required=True, type=str)
    parser.add_argument("--out", help="Output file name. Default: stdout", required=False, type=str)
    parser.add_argument("--region", help="Region chr:start-end", required=False, type=str)
    args = parser.parse_args()

    # Prepare output
    if args.out: outf = open(args.out,"w")
    else: outf = sys.stdout

    vcf = VCF(args.vcf) #vcf.Reader(open(args.vcf, "rb"))
    if args.region:
        vcf = vcf(args.region) #reader.fetch(args.region)

    for record in vcf:
        genotypes = np.array(record.genotypes)[:,0:2]

        y_val = np.apply_along_axis(gp_position_np, 1, genotypes)
        x_val = range(y_val.shape[0])

        output = " ".join([str(record.CHROM), str(record.POS), str(record.ID), str(np.mean(record.format('GP')[x_val, y_val]))])
        PrintLine(output, outf)

if __name__ == "__main__":
    main()
