#!/usr/bin/env python
"""

Author: Shubham Saini
shubhamsaini@eng.ucsd.edu

"""


import argparse
import sys

from cyvcf2 import VCF, Writer
from itertools import combinations_with_replacement
import numpy as np
from collections import Counter
import pandas as pd
import seaborn as sns

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

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--vcf", help="Input VCF file", required=True, type=str)
    parser.add_argument("--chrom", help="Chromosome", required=True, type=int)
    parser.add_argument("--position", help="Position", required=True, type=int)
    args = parser.parse_args()

    reader = VCF(args.vcf)
    vcf_samples = reader.samples
    reader = reader(f'{args.chrom}:{args.position}-{args.position}')
    iqr_outliers = True #
    iqr_outliers_min_samples = 100 #
    rmrare = 0.05

    for record in reader:
        gtdata = {}
        exclude_samples = []
        alleles = [record.REF]+record.ALT
        genotypes = np.array(record.genotypes)[:,0:2]

        acounts = Counter(genotypes.flatten())
        afreqs = [acounts[i]/genotypes.flatten().shape[0] for i in range(len(alleles))]
        aaf = sum(afreqs[1:])
        aaf = min([aaf, 1-aaf])

        alleles_lengths = [len(record.REF)]+[len(i) for i in record.ALT]
        comb = list(combinations_with_replacement(range(len(alleles_lengths)), 2))
        idx = [gp_position(genotype[0], genotype[1]) for genotype in comb]
        geno_sum_lengths = [(alleles_lengths[genotype[0]] + alleles_lengths[genotype[1]]) for genotype in comb]
        geno_sum_lengths_dict = dict(zip(idx, geno_sum_lengths))
        geno_sum_lengths_sorted = [geno_sum_lengths_dict[i] for i in sorted(geno_sum_lengths_dict)]
        #geno_sum_lengths_sorted = [geno_sum_lengths[i] for i in idx]

        gp_sum = np.dot(geno_sum_lengths_sorted, record.format('GP').T)
        gtdata = dict(zip(vcf_samples, gp_sum))

        if iqr_outliers:
            q75, q25 = np.percentile(gp_sum, [75 ,25])
            iqr = q75 - q25
            low_threshold = q25 - (1.5*iqr)
            high_threshold = q75 + (1.5*iqr)

            if iqr_outliers_min_samples > 0:
                alleles_idx = list(range(len(alleles_lengths)))
                genotypes_lenghts = replace_allele_idx(genotypes, alleles_idx, alleles_lengths)
                genotypes_lenghts = np.sum(genotypes_lenghts, axis=1)
                gt_counts = Counter(genotypes_lenghts)
                valid_gt_counts = [i for i in gt_counts if gt_counts[i] > iqr_outliers_min_samples]
                if len(valid_gt_counts) > 0:
                    min_gt_counts = min(valid_gt_counts)
                    max_gt_counts = max(valid_gt_counts)
                    low_threshold = min(low_threshold, min_gt_counts)
                    high_threshold = max(high_threshold, max_gt_counts)

            exclude_samples = [sample for sample in gtdata if (gtdata[sample]>=high_threshold or gtdata[sample]<=low_threshold)]

        dosage_data = pd.DataFrame.from_dict(gtdata, orient='index').reset_index()
        dosage_data.columns = ['sample','dosage']
        dosage_data['exclude'] = 0
        dosage_data.loc[dosage_data['sample'].isin(exclude_samples),'exclude'] = 1
        print("Excluded", len(exclude_samples), "samples")

        gtdata = {}
        exclude_samples = []
        alleles = [record.REF]+record.ALT
        genotypes = np.array(record.genotypes)[:,0:2]

        acounts = Counter(genotypes.flatten())
        afreqs = [acounts[i]/genotypes.flatten().shape[0] for i in range(len(alleles))]
        aaf = sum(afreqs[1:])
        aaf = min([aaf, 1-aaf])

        for sample in range(len(vcf_samples)):
            f1, f2 = [afreqs[int(item)] for item in genotypes[sample]]
            if f1 < rmrare or f2 < rmrare:
                exclude_samples.append(vcf_samples[sample])
                gtdata[vcf_samples[sample]] = sum([len(record.REF) for item in genotypes[sample]])
            else:
                gtdata[vcf_samples[sample]] = sum([len(alleles[int(item)]) for item in genotypes[sample]])

        genotype_data = pd.DataFrame.from_dict(gtdata, orient='index').reset_index()
        genotype_data.columns = ['sample','genotype']
        genotype_data.loc[genotype_data['sample'].isin(exclude_samples),'genotype'] = np.nan

    genotype_data = genotype_data.merge(dosage_data)

    xvals = genotype_data['genotype'].values
    yvals = genotype_data['dosage'].values
    combos = list(zip(xvals, yvals))
    weight_counter = Counter(combos)
    genotype_data['freq'] = genotype_data.apply(lambda x: weight_counter[(x['genotype'], x['dosage'])], axis="columns")

    sns_plot = sns.scatterplot(data=genotype_data, x="genotype", y="dosage", hue="freq", size="freq", sizes=(20, 200))
    plt.axhline(low_threshold, color='red', linestyle='dashed')
    plt.axhline(high_threshold, color='red', linestyle='dashed')
    fig = sns_plot.get_figure()
    fig.savefig(f'chr{args.chrom}_{args.position}_dosage_vs_gt.png')

if __name__ == "__main__":
    main()
