#!/usr/bin/env python
"""

Author: Shubham Saini
shubhamsaini@eng.ucsd.edu

"""


import argparse
import sys
import numpy as np
from cyvcf2 import VCF
from scipy.stats import binom_test
from collections import Counter

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

    vcf = VCF(args.vcf)
    if args.region:
        vcf = vcf(args.region)

    output = f'CHROM\tPOS\tSNP\tALLELE\tAFREQ'
    PrintLine(output, outf)

    for record in vcf:
        alleles = [record.REF]+record.ALT
        genotypes = np.array(record.genotypes)[:,0:2]
        acounts = Counter(genotypes.flatten())
        afreqs = [acounts[i]/genotypes.flatten().shape[0] for i in range(len(alleles))]

        for i in range(len(alleles)):
            output = f'{record.CHROM}\t{record.POS}\t{record.ID}\t{alleles[i]}\t{afreqs[i]}'
            PrintLine(output, outf)

if __name__ == "__main__":
    main()
