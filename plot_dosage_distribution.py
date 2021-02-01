#!/usr/bin/env python
"""

Author: Shubham Saini
shubhamsaini@eng.ucsd.edu

"""


import argparse
import sys
import numpy as np
from cyvcf2 import VCF
from itertools import combinations_with_replacement
import matplotlib.pyplot as plt

def gp_position(j, k):
    p1 = (k*(k+1)/2)+j
    p2 = (j*(j+1)/2)+k
    return int(max(p1, p2))

def PrintLine(text, f):
    f.write(text+"\n")
    f.flush()

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--vcf", help="Argument 1 description", required=True, type=str)
    parser.add_argument("--chrom", help="Region chr:start-end", required=True, type=int)
    parser.add_argument("--pos", help="Region chr:start-end", required=True, type=int)
    args = parser.parse_args()

    vcf = VCF(args.vcf)
    region = '%d:%d-%d'%(args.chrom, args.pos, args.pos)

    for record in vcf(region):
        genotypes = np.array(record.genotypes)[:,0:2]

        alleles_lengths = [len(record.REF)]+[len(i) for i in record.ALT]
        comb = list(combinations_with_replacement(range(len(alleles_lengths)), 2))
        idx = [gp_position(genotype[0], genotype[1]) for genotype in comb]
        geno_sum_lengths = [(alleles_lengths[genotype[0]] + alleles_lengths[genotype[1]]) for genotype in comb]
        geno_sum_lengths_sorted = [geno_sum_lengths[i] for i in idx]
        gp_sum = np.sort(np.dot(geno_sum_lengths_sorted, record.format('GP').T))

        #plt.scatter(range(len(gp_sum)), gp_sum, marker='.')
        #plt.xlabel("Sample")
        #plt.ylabel("Dosage");

        plt.hist(gp_sum, bins='auto')
        # Mean +/- 3 SD
        #plt.axvline(np.mean(gp_sum)+(3*np.std(gp_sum)), color='red', linestyle='dashed')
        #plt.axvline(np.mean(gp_sum)-(3*np.std(gp_sum)), color='red', linestyle='dashed')
        #plt.xlabel("Dosage")
        #min_ylim, max_ylim = plt.ylim()
        #plt.text((np.mean(gp_sum)+(3*np.std(gp_sum)))+0.01, max_ylim*0.95, 'Mean + 3 S.D.')
        #plt.text((np.mean(gp_sum)-(3*np.std(gp_sum)))+0.01, max_ylim*0.95, 'Mean - 3 S.D.')

        # IQR filtering
        # Q1 - 1.5IQR; Q3 + 1.5IQR
        q75, q25 = np.percentile(gp_sum, [75 ,25])
        iqr = q75 - q25
        plt.axvline(q75 + (1.5*iqr), color='red', linestyle='dashed')
        plt.axvline(q25 - (1.5*iqr), color='red', linestyle='dashed')
        plt.xlabel("Dosage")
        min_ylim, max_ylim = plt.ylim()
        plt.text((q75 + (1.5*iqr))+0.01, max_ylim*0.95, 'Q1 - 1.5IQR')
        plt.text((q25 - (1.5*iqr))+0.01, max_ylim*0.95, 'Q3 + 1.5IQR')

        plt.title('chr%d:%d'%(args.chrom, args.pos))
        plt.savefig('chr%d.%d.dosage.png'%(args.chrom, args.pos))

if __name__ == "__main__":
    main()
