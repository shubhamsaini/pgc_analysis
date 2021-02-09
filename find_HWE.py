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

def replace_allele_idx(numbers, problem_numbers, alternative_numbers):
    numbers = np.asarray(numbers)
    problem_numbers = np.asarray(problem_numbers)
    alternative_numbers = np.asarray(alternative_numbers)
    n_min, n_max = numbers.min(), numbers.max()
    replacer = np.arange(n_min, n_max + 1)
    # Mask replacements out of range
    mask = (problem_numbers >= n_min) & (problem_numbers <= n_max)
    replacer[problem_numbers[mask] - n_min] = alternative_numbers[mask]
    numbers = replacer[numbers - n_min]
    return numbers

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
        alleles = [record.REF]+record.ALT
        alleles_lengths = [len(record.REF)]+[len(i) for i in record.ALT]
        alleles_idx = list(range(len(alleles_lengths)))
        genotypes = np.array(record.genotypes)[:,0:2]
        genotypes_lenghts = replace_allele_idx(genotypes, alleles_idx, alleles_lengths)

        acounts = Counter(genotypes_lenghts.flatten())
        afreqs = [i/sum(acounts.values()) for i in acounts.values()]

        exp_hom_frac = sum([val**2 for val in afreqs])
        num_hom = sum(genotypes_lenghts[:,0] == genotypes_lenghts[:,1])
        total_samples = genotypes_lenghts.shape[0]

        hwe_p = binom_test(num_hom, n=total_samples, p=exp_hom_frac)

        output = " ".join([str(record.CHROM), str(record.POS), str(record.ID), str(hwe_p)])
        PrintLine(output, outf)

if __name__ == "__main__":
    main()
