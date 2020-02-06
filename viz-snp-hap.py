#!/usr/bin/env python
"""

Author: Shubham Saini
shubhamsaini@eng.ucsd.edu

Example:
./viz-snp-hap.py --vcf /project/gymrek/chr10/scz_munc_eur-qc.bgs.chr10.imputed.vcf.gz --str-id STR_187806 --pos 104639652 --chrom 10 --samples-cases cases.txt --samples-controls controls.txt --out-dir test_viz_out/

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
    parser.add_argument("--pos", help="STR Genomic Coordinate", required=True, type=int)
    parser.add_argument("--chrom", help="STR Chromosome", required=True, type=int)
    parser.add_argument("--window", help="Haplotype window size", required=False, type=int, default=5000)
    parser.add_argument("--out-dir", help="Figure output directory", required=False, type=str, default = ".")
    parser.add_argument("--samples-cases", help="Cases file", required=True, type=str)
    parser.add_argument("--samples-controls", help="Controls file", required=True, type=str)
    args = parser.parse_args()

    str_id = args.str_id #"STR_187806"
    START = args.pos #104639652
    CHROM= args.chrom #10
    WINDOW = args.window

    vcf_files = []
    with open(args.vcf) as f:
        for line in f:
            vcf_files.append(line.strip())

    snp_info = {}
    str_info = {}
    str_gt = {}
    snp_gt = {}
    all_samples = []

    for vcf_file_loc in vcf_files:
        if args.samples_cases:
            with open(args.samples_cases,"r") as samples_file:
                samples_cases = [line.rstrip('\n') for line in samples_file]
            with open(args.samples_controls,"r") as samples_file:
                samples_controls = [line.rstrip('\n') for line in samples_file]
            samples = samples_cases + samples_controls
            vcf = VCF(vcf_file_loc, samples=samples)
        else:
            vcf = VCF(vcf_file_loc)

        all_samples = all_samples + vcf.samples
        for v in vcf("%s:%s-%s"%(CHROM, START-WINDOW, START+WINDOW)):
            if v.ID == str_id:
                if v.ID not in str_info.keys():
                    str_info[v.ID] = [v.POS, v.REF, ",".join(v.ALT)]
                refLen = len(v.REF)
                gt_bases = v.gt_bases
                str_dosages = []
                for gt in gt_bases:
                    if ('.' in gt):
                        str_dosages = str_dosages + [np.nan, np.nan]
                    else:
                        gt_bases_split = re.split('/|\|',gt)
                        str_dosages = str_dosages + [(len(i) - refLen) for i in gt_bases_split]
                if v.ID not in str_gt.keys():
                    str_gt[v.ID] = str_dosages
                else:
                    str_gt[v.ID] = str_gt[v.ID] + str_dosages
            else:
                refLen = len(v.REF)
                if refLen > 1:
                    continue
                if len(v.ALT) > 1:
                    continue

                if v.ID not in snp_info.keys():
                    snp_info[v.ID] = [v.POS, v.REF, ",".join(v.ALT)]
                snp_dosages = []
                for x in v.genotypes:
                    if ('.' in x):
                        snp_dosages = snp_dosages + [np.nan, np.nan]
                    else:
                        snp_dosages = snp_dosages + [x[0], x[1]]
                if v.ID not in snp_gt.keys():
                    snp_gt[v.ID] = snp_dosages
                else:
                    snp_gt[v.ID] = snp_gt[v.ID] + snp_dosages


    str_info = pd.DataFrame.from_dict(str_info, orient="index", columns=['pos','ref','alt'])
    snp_info = pd.DataFrame.from_dict(snp_info, orient="index", columns=['pos','ref','alt'])
    colnames = []
    for i in all_samples:
        colnames.extend([i+"_0", i+"_1"])
    snp_gt = pd.DataFrame.from_dict(snp_gt, orient="index", columns=colnames)
    str_gt = pd.DataFrame.from_dict(str_gt, orient="index", columns=colnames)
    snp_gt.dropna(inplace=True)
    str_gt.dropna(inplace=True)
    str_info = str_info.merge(str_gt, left_index=True, right_index=True)
    snp_info = snp_info.merge(snp_gt, left_index=True, right_index=True)
    haplotypes_full = snp_info.append(str_info).reset_index()
    haplotypes_full.rename(columns={"index": "id"}, inplace=True)
    print(str_info.shape)
    print(snp_info.shape)

    snp = snp_gt.values.T
    str_dosages = str_gt.values.tolist()[0]

    elastic_net= ElasticNet(alpha= 0.1, l1_ratio= 0.5) # l1_ratio is the mix r
    elastic_net.fit(snp,str_dosages)
    print("Model accuracy: ", elastic_net.score(snp,str_dosages))

    HAP_LENGTH = 20
    locus = str_id.replace("/","_")
    RSIDS=[snp_gt.index.tolist()[i] for i in np.argsort(np.absolute(elastic_net.coef_))[-HAP_LENGTH:]] + [locus]
    print("Top SNPS:")
    print(RSIDS[:-1])

    samples_cases = list(set(samples_cases).intersection(all_samples))
    samples_controls = list(set(samples_controls).intersection(all_samples))


    cases_colnames = []
    for i in samples_cases:
        cases_colnames.extend([i+"_0", i+"_1"])
    cases_colnames = ["id","pos","ref","alt"] + cases_colnames


    controls_colnames = []
    for i in samples_controls:
        controls_colnames.extend([i+"_0", i+"_1"])
    controls_colnames = ["id","pos","ref","alt"] + controls_colnames

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
