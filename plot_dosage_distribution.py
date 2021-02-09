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
import pandas as pd
from collections import Counter

def gp_position(j, k):
    p1 = (k*(k+1)/2)+j
    p2 = (j*(j+1)/2)+k
    return int(max(p1, p2))

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
    parser.add_argument("--vcf", help="VCF File", required=True, type=str)
    parser.add_argument("--chrom", help="Region chr:start-end", required=True, type=int)
    parser.add_argument("--pos", help="Region chr:start-end", required=True, type=int)
    parser.add_argument("--min-gt-counts", help="Min GT counts", required=False, type=int, default=0)
    parser.add_argument("--fam", help="FAM file", required=False, type=str)
    parser.add_argument("--out-prefix", help="Output Prefix", required=False, type=str)
    args = parser.parse_args()

    vcf = VCF(args.vcf)
    region = '%d:%d-%d'%(args.chrom, args.pos, args.pos)

    if args.fam:
        fam_data = pd.read_csv(args.fam, names = ['FID','IID','FATHER','MOTHER','SEX','PHENO'], delim_whitespace=True)
        fam_data['SAMPLE'] = fam_data.apply(lambda x: x["FID"]+"_"+x["IID"], 1)
        fam_data = fam_data[['SAMPLE','PHENO']]

        sorter = vcf.samples
        sorterIndex = dict(zip(sorter, range(len(sorter))))
        fam_data['sample_rank'] = fam_data['SAMPLE'].map(sorterIndex)
        fam_data.sort_values(['sample_rank'], ascending = True, inplace = True)
        fam_data.drop('sample_rank', 1, inplace = True)

        fam = list(fam_data['PHENO'].values)
    else:
        fam = [0]*len(vcf.samples)

    for record in vcf(region):
        genotypes = np.array(record.genotypes)[:,0:2]

        alleles_lengths = [len(record.REF)]+[len(i) for i in record.ALT]
        comb = list(combinations_with_replacement(range(len(alleles_lengths)), 2))
        idx = [gp_position(genotype[0], genotype[1]) for genotype in comb]
        geno_sum_lengths = [(alleles_lengths[genotype[0]] + alleles_lengths[genotype[1]]) for genotype in comb]
        geno_sum_lengths_sorted = [geno_sum_lengths[i] for i in idx]
        gp_sum = np.sort(np.dot(geno_sum_lengths_sorted, record.format('GP').T))

        q75, q25 = np.percentile(gp_sum, [75 ,25])
        iqr = q75 - q25
        low_threshold = q25 - (1.5*iqr)
        high_threshold = q75 + (1.5*iqr)

        if args.min_gt_counts > 0:
            alleles_idx = list(range(len(alleles_lengths)))
            genotypes_lenghts = replace_allele_idx(genotypes, alleles_idx, alleles_lengths)
            genotypes_lenghts = np.sum(genotypes_lenghts, axis=1)
            gt_counts = Counter(genotypes_lenghts)
            valid_gt_counts = [i for i in gt_counts if gt_counts[i] > args.min_gt_counts]
            min_gt_counts = min(valid_gt_counts)
            max_gt_counts = max(valid_gt_counts)
            low_threshold = min(low_threshold, min_gt_counts)
            high_threshold = max(high_threshold, max_gt_counts)

        df = pd.DataFrame({"x" : gp_sum,"class" : fam})

        _, edges = np.histogram(df["x"], bins='auto')
        histdata = []; labels=[]
        for n, group in df.groupby("class"):
            histdata.append(np.histogram(group["x"], bins=edges)[0])
            labels.append(n)

        hist = np.array(histdata)
        histcum = np.cumsum(hist,axis=0)

        plt.bar(edges[:-1],hist[0,:], width=np.diff(edges)[0],
                    label=labels[0], align="edge")

        for i in range(1,len(hist)):
            plt.bar(edges[:-1],hist[i,:], width=np.diff(edges)[0],
                    bottom=histcum[i-1,:],label=labels[i], align="edge")


        plt.axvline(low_threshold, color='red', linestyle='dashed')
        plt.axvline(high_threshold, color='red', linestyle='dashed')
        plt.xlabel("Dosage")
        min_ylim, max_ylim = plt.ylim()
        plt.text(high_threshold+0.01, max_ylim*0.95, 'Q3 + 1.5IQR')
        plt.text(low_threshold+0.01, max_ylim*0.95, 'Q1 - 1.5IQR')

        plt.legend(title="class")
        plt.title('chr%d:%d'%(args.chrom, args.pos))
        if args.out_prefix:
            plt.savefig('%s.chr%d.%d.dosage.pheno.png'%(args.out_prefix,args.chrom, args.pos))
        else:
            plt.savefig('chr%d.%d.dosage.pheno.png'%(args.chrom, args.pos))

if __name__ == "__main__":
    main()
