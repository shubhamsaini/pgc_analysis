#!/usr/bin/env python
"""

Author: Shubham Saini
shubhamsaini@eng.ucsd.edu

Example:
./viz-snp-hap.py --vcf /project/gymrek/chr10/scz_munc_eur-qc.bgs.chr10.imputed.vcf.gz --str-id STR_187806 --pos 10:104639652-104639652 --samples-cases cases.txt --samples-controls controls.txt --out-dir hap_figs/

"""

# Allow us to edit fonts in Illustrator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os
import scipy.stats
import seaborn as sns
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib import gridspec


from cyvcf2 import VCF
import sys
import re
import numpy as np
from sklearn.linear_model import ElasticNet

import argparse
import sys

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def PlotHapmaptrix(hapmatrix, allele, allsnps, fname):
    box_w =  1.0/len(allsnps)
    box_h = box_w
    hap_height = hapmatrix.shape[0]*0.0025*4
    legend_height = 0
    fig = plt.figure()
    fig.set_size_inches(3, hap_height + legend_height)
    #gs = gridspec.GridSpec(2, 1, height_ratios=[hap_height, legend_height])
    ax = fig.add_subplot(111)
    # Plot SNPs
    imx = ax.imshow(hapmatrix, vmin=0, vmax=1, cmap=plt.cm.Greys.from_list("snp", ["lightgray","black"]),
              aspect="auto", extent=(0, hapmatrix.shape[1], box_h, hapmatrix.shape[0]-box_h))
    ax.set_yticks([]);
    ax.set_yticklabels([]);
    ax.set_xticks([]);
    ax.set_xticklabels([]);
    ax.set_title("STR allele %s"%allele)
    fig.savefig(fname)

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--vcf", help="VCF Input", required=True, type=str)
    parser.add_argument("--str-id", help="STR ID", required=True, type=str)
    parser.add_argument("--pos", help="STR POSITION (chr:start-end)", required=True, type=str)
    parser.add_argument("--out-dir", help="Figure output directory", required=True, type=str)
    parser.add_argument("--samples-cases", help="Cases file", required=False, type=str)
    parser.add_argument("--samples-controls", help="Controls file", required=False, type=str)
    args = parser.parse_args()

    hap_file = args.vcf #"/storage/s1saini/hipstr_allfilters/str_snp/chr10.str.snp.feb18.vcf.gz"
    str_id = args.str_id #"STR_187806"
    POS = args.pos #"10:104639652-104639652"

    haplotype_df = []

    if args.samples_cases:
        with open(args.samples_cases,"r") as samples_file:
            samples_cases = [line.rstrip('\n') for line in samples_file]
        with open(args.samples_controls,"r") as samples_file:
            samples_controls = [line.rstrip('\n') for line in samples_file]
        samples = samples_cases + samples_controls
        vcf = VCF(hap_file, samples=samples)
    else:
        vcf = VCF(hap_file)
    gt_bases_len = []
    for v in vcf(POS):
        if v.ID != str_id:
            continue
        refLen = len(v.REF)
        if refLen == 1:
            break

        curr_hap = [v.ID, v.POS, v.REF, ",".join(v.ALT)]
        CHROM = int(v.CHROM)
        START = int(v.POS)
        gt_bases = v.gt_bases
        for gt in gt_bases:
            if ('.' in gt):
                gt_bases_len = gt_bases_len + [np.nan, np.nan]
                curr_hap = curr_hap + [np.nan, np.nan]
            else:
                gt_bases_split = re.split('/|\|',gt)
                gt_bases_len = gt_bases_len + [(len(i) - refLen) for i in gt_bases_split]
                curr_hap = curr_hap + [(len(i) - refLen) for i in gt_bases_split]


    str_dosages = gt_bases_len

    numhaps = len(vcf.samples)*2
    colnames = ["id","pos","ref","alt"] + ["hap_%s"%i for i in range(numhaps)]
    haplotype_df.append(curr_hap)
    print(len(haplotype_df))

    WINDOW = 100000

    if args.samples_cases:
        with open(args.samples_cases,"r") as samples_file:
            samples_cases = [line.rstrip('\n') for line in samples_file]
        with open(args.samples_controls,"r") as samples_file:
            samples_controls = [line.rstrip('\n') for line in samples_file]
        samples = samples_cases + samples_controls
        vcf = VCF(hap_file, samples=samples)
    else:
        vcf = VCF(hap_file)
    snp_gt = []
    snp_id = []
    for v in vcf("%s:%s-%s"%(CHROM, START-WINDOW, START+WINDOW)):
        refLen = len(v.REF)
        if refLen > 1:
            continue
        if len(v.ALT) > 1:
            continue

        curr_hap = [v.ID, v.POS, v.REF, ",".join(v.ALT)]
        snp_id.append(str(v.ID))
        gt_bases_len = []
        for x in v.genotypes:
            if ('.' in x):
                gt_bases_len = gt_bases_len + [np.nan, np.nan]
                curr_hap = curr_hap + [np.nan, np.nan]
            else:
                gt_bases_len = gt_bases_len + [x[0], x[1]]
                curr_hap = curr_hap + [x[0], x[1]]

        haplotype_df.append(curr_hap)
        snp_gt.append(gt_bases_len)

    snp = np.array(snp_gt).T
    haplotypes_full = pd.DataFrame(haplotype_df, columns=colnames)

    print(snp.shape)
    print(len(str_dosages))

    print(haplotypes_full.shape)

    elastic_net= ElasticNet(alpha= 0.1, l1_ratio= 0.5) # l1_ratio is the mix r

    elastic_net.fit(snp,str_dosages)
    print(elastic_net.score(snp,str_dosages))

    HAP_LENGTH = 20
    locus = str_id.replace("/","_")
    RSIDS=[snp_id[i] for i in np.argsort(np.absolute(elastic_net.coef_))[-HAP_LENGTH:]] + [locus]
    print("Top SNPS:")
    print(RSIDS)

    samples_cases = list(set(samples_cases).intersection(vcf.samples))
    samples_controls = list(set(samples_controls).intersection(vcf.samples))
    num_cases = len(samples_cases)
    num_controls = len(samples_controls)

    cases_colnames = ["id","pos","ref","alt"] + ["hap_%s"%i for i in range(num_cases*2)]
    controls_colnames = ["id","pos","ref","alt"] + ["hap_%s"%i for i in range(num_cases*2+1, numhaps)]

    ####
    ### Do for cases first
    ####
    print("Viz cases")
    haplotypes = haplotypes_full[cases_colnames]
    haplotypes = haplotypes[haplotypes['id'].isin(RSIDS)]
    haplotypes["vartype"] = haplotypes.apply(lambda x: ["SNP","STR"][int(len(x["ref"])>1)], 1)
    haplotypes.index = [a+b for a, b in zip(haplotypes["vartype"].values, [str(i) for i in range(haplotypes.shape[0])] )]

    # Annotate STR lengths
    ref = haplotypes[haplotypes["vartype"]=="STR"]["ref"].values[0]
    alt = haplotypes[haplotypes["vartype"]=="STR"]["alt"].values[0].split(",")
    str_allele_lengths = [len(ref)] + [len(item) for item in alt]
    str_allele_lengths = [item-len(ref) for item in str_allele_lengths]

    hapcols = cases_colnames[4:]
    haplotype_filt = haplotypes[hapcols].transpose()

    allsnps = [item for item in haplotype_filt.columns if "SNP" in item]
    str_col = list(set(haplotype_filt.columns) - set(allsnps))[0]
    haplotype_filt = haplotype_filt.sort_values(by=str_col)

    for allele in sorted(list(set(str_allele_lengths)))+["max"]:
        if allele == "max":
            hapmatrix = np.matrix(haplotype_filt[haplotype_filt[str_col]>=15][allsnps])
        else: hapmatrix = np.matrix(haplotype_filt[haplotype_filt[str_col]==allele][allsnps])
        if hapmatrix.shape[0]>= 10:
            sys.stderr.write("%s:%s\n"%(allele, hapmatrix.shape))
            fname = os.path.join("./", "%s/Cases_Haplotypes_%s.pdf"%(args.out_dir,allele))
            PlotHapmaptrix(hapmatrix, allele, allsnps, fname)


    ###
    ### Do for controls now
    ####
    print("Viz controls")
    haplotypes = haplotypes_full[controls_colnames]
    haplotypes = haplotypes[haplotypes['id'].isin(RSIDS)]
    haplotypes["vartype"] = haplotypes.apply(lambda x: ["SNP","STR"][int(len(x["ref"])>1)], 1)
    haplotypes.index = [a+b for a, b in zip(haplotypes["vartype"].values, [str(i) for i in range(haplotypes.shape[0])] )]

    # Annotate STR lengths
    ref = haplotypes[haplotypes["vartype"]=="STR"]["ref"].values[0]
    alt = haplotypes[haplotypes["vartype"]=="STR"]["alt"].values[0].split(",")
    str_allele_lengths = [len(ref)] + [len(item) for item in alt]
    str_allele_lengths = [item-len(ref) for item in str_allele_lengths]

    hapcols = controls_colnames[4:]
    haplotype_filt = haplotypes[hapcols].transpose()

    allsnps = [item for item in haplotype_filt.columns if "SNP" in item]
    str_col = list(set(haplotype_filt.columns) - set(allsnps))[0]
    haplotype_filt = haplotype_filt.sort_values(by=str_col)

    for allele in sorted(list(set(str_allele_lengths)))+["max"]:
        if allele == "max":
            hapmatrix = np.matrix(haplotype_filt[haplotype_filt[str_col]>=15][allsnps])
        else: hapmatrix = np.matrix(haplotype_filt[haplotype_filt[str_col]==allele][allsnps])
        if hapmatrix.shape[0]>= 10:
            sys.stderr.write("%s:%s\n"%(allele, hapmatrix.shape))
            fname = os.path.join("./", "%s/Controls_Haplotypes_%s.pdf"%(args.out_dir,allele))
            PlotHapmaptrix(hapmatrix, allele, allsnps, fname)

if __name__ == "__main__":
    main()
