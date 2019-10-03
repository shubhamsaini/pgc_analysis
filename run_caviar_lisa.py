#!/usr/bin/env python
"""

Author: Shubham Saini
shubhamsaini@eng.ucsd.edu


Run CAVIAR on SNP-STR GWAS results

Example usage:
./run_caviar_lisa.py \
--str-vcf /nfs/gymrek/chr22/strs/pgc.strs.imputed.chr22.vcf.gz \
--snp-vcf /home/gymrek/ssaini/caviar/chr22.files.list \
--str-assoc /home/gymrek/ssaini/gwas/assoc_results/str_assoc/pgc.chr22.assoc \
--snp-assoc /home/gymrek/ssaini/gwas/assoc_results/snp_assoc/chr22.snps.metal.txt \
--snp-rsid /home/gymrek/ssaini/gwas/assoc_results/snp_assoc/chr22.pos.rsid.txt \
--chrom 22 --pos 39975317 \
--window-kb 10

"""
import matplotlib
matplotlib.use('Agg')

import argparse
import sys

from cyvcf2 import VCF
import numpy as np
import pandas as pd
from scipy.stats import norm
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
import re

from subprocess import Popen, PIPE, STDOUT

try:
    from subprocess import DEVNULL # py3k
except ImportError:
    import os
    DEVNULL = open(os.devnull, 'wb')

def PrintLine(text, f):
    f.write(text+"\n")
    f.flush()


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--str-vcf", help="STR VCF file for LD calculation", required=True, type=str)
    parser.add_argument("--snp-vcf", help="SNP VCF files list for LD calculation", required=True, type=str)
    parser.add_argument("--str-assoc", help="STR GWAS file from plinkSTR", required=True, type=str)
    parser.add_argument("--snp-assoc", help="SNP GWAS file from Metal", required=True, type=str)
    parser.add_argument("--snp-rsid", help="SNP GWAS file from Metal", required=True, type=str)
    parser.add_argument("--chrom", help="CHROM", required=True, type=str)
    parser.add_argument("--pos", help="POSITION", required=True, type=str)
    parser.add_argument("--window-kb", help="Window in KB", required=True, type=str)
    parser.add_argument("--maxp", help="Maximum p-value", required=False, type=str)

    parser.add_argument("--out", help="Output file name. Default: stdout", required=False, type=str)
    args = parser.parse_args()

    # Prepare output
    if args.out: outf = open(args.out,"w")
    else: outf = sys.stdout

    snp_vcf_files = []
    with open(args.snp_vcf) as f:
        for line in f:
            snp_vcf_files.append(line.strip())

    str_vcf_file = args.str_vcf

    CHROM = int(args.chrom)
    POS = int(args.pos)
    window=int(args.window_kb)*1000
    START = POS-window
    END = POS+window

    ### CAVIAR parameters
    seed = np.random.randint(100000)
    ldfile = 'caviar_files/corr_mat_%d.txt'%(seed)
    zfile = 'caviar_files/zscore_%d.csv'%(seed)
    numcausal = "2"
    outfile = "caviar_files/caviar.output_%d_%d_%d"%(CHROM, POS, seed)


    PrintLine('Seed %d %d %d %d'%(seed, CHROM, START, END), outf)
    PrintLine("Loading STR GWAS results", outf)

    str_assoc_result = pd.read_csv(args.str_assoc, delim_whitespace=True)
    str_assoc_result = str_assoc_result[[str_assoc_result.columns.tolist()[i] for i in [0,1,3,2,6,4]]]
    str_assoc_result.columns = ["CHR", "BP","P","ID","CI95","OR"]
    str_assoc_result = str_assoc_result[(str_assoc_result['BP']>=START) & (str_assoc_result['BP']<=END)]
    str_assoc_result = str_assoc_result[str_assoc_result.P != "NAN"]
    str_assoc_result['vtype'] = "STR"
    str_assoc_result['P'] = str_assoc_result.P.astype('float')

    if args.maxp:
        str_assoc_result = str_assoc_result[str_assoc_result['P']<=float(args.maxp)]
    str_assoc_result = str_assoc_result.sort_values(by=['BP'])
    print(str_assoc_result.shape)

    PrintLine("Loading SNP GWAS results", outf)

    snp_assoc_result = pd.read_csv(args.snp_assoc, delim_whitespace=True)
    snp_assoc_result = snp_assoc_result[snp_assoc_result['TEST']=="ADD"]
    snp_assoc_result = snp_assoc_result[["POS","P","ID","Z_STAT"]]
    snp_assoc_result.columns = ["BP","P","ID","Zscore"]
    snp_assoc_result['BP'] = snp_assoc_result.BP.astype('int')
    snp_assoc_result = snp_assoc_result[(snp_assoc_result['BP']>=START) & (snp_assoc_result['BP']<=END)]
    snp_assoc_result['vtype'] = "SNP"
    snp_assoc_result['CHR'] = CHROM
    snp_assoc_result['P'] = snp_assoc_result.P.astype('float')
    if args.maxp:
        snp_assoc_result = snp_assoc_result[snp_assoc_result['P']<=float(args.maxp)]
    snp_assoc_result = snp_assoc_result.sort_values(by=['BP'])
    print(snp_assoc_result.shape)

    gwas_rsid = list(snp_assoc_result['ID'].values) + list(str_assoc_result['ID'].values)


    ### STR samples
    vcf = VCF(str_vcf_file, gts012=True)
    str_samples = vcf.samples


    region = '%d:%d-%d'%(CHROM,START,END)

    for vcf_file_loc in snp_vcf_files:
        cmd = "bcftools query -r "+region+" -f '%ID\n' "+vcf_file_loc+" >> caviar_files/snp_files_rsid_" + str(seed) + ".txt"
        p = Popen(cmd, shell=True, stdout=DEVNULL, stderr=DEVNULL)
        output = p.communicate()[0]
        if p.returncode != 0:
            print("Failed")

    snp_all_rsid = set()
    with open("caviar_files/snp_files_rsid_%d.txt"%(seed)) as f:
        for line in f:
            if line.strip() in gwas_rsid:
                snp_all_rsid.add(line.strip())


    PrintLine("Loading SNP genotypes", outf)
    snps_rsid = set()
    str_rsid = set()
    samples = []

    gt_array = defaultdict(list)
    for rsid in snp_all_rsid:
        if rsid in gwas_rsid:
            gt_array[rsid] = []

    for vcf_file_loc in snp_vcf_files:
        vcf_rsid = set()
        vcf = VCF(vcf_file_loc, gts012=True, samples = str_samples)
        for v in vcf('%d:%d-%d'%(CHROM,START,END)):
            vcf_rsid.add(str(v.ID))
            if len(v.REF) == 1 and str(v.ID) in gwas_rsid:
                gt = np.array(v.gt_types, dtype="float")
                gt[gt==3] = np.nan
                gt_array[str(v.ID)] = gt_array[str(v.ID)] + list(gt)
                snps_rsid.add(str(v.ID))

        samples = samples + list(vcf.samples)
        not_found_rsid = snp_all_rsid.difference(vcf_rsid)
        for rsid in not_found_rsid:
            gt_array[rsid] = gt_array[rsid] + [np.nan]*len(vcf.samples)

    for i in gt_array.keys():
        if len(gt_array[i]) == 0:
            del gt_array[i]

    #print([len(gt_array[i]) for i in gt_array.keys()])
    gt_array = pd.DataFrame(gt_array, dtype="float")
    print(gt_array.shape)

    PrintLine("Loading STR genotypes",outf)
    gt_array2 = defaultdict(list)

    vcf = VCF(str_vcf_file, gts012=True, samples=samples)
    for v in vcf('%d:%d-%d'%(CHROM,START,END)):
        str_gt = []
        if str(v.ID) in gwas_rsid:
            for gt in v.gt_bases:
                str_gt.append(np.sum([len(i)-len(v.REF) for i in re.split('/|\|',gt)]))
            gt_array2[str(v.ID)] = str_gt
            str_rsid.add(str(v.ID))

    gt_array2 = pd.DataFrame(gt_array2, dtype="float")
    print(gt_array2.shape)

    PrintLine("Writing LD file", outf)
    full_geno_matrix = pd.concat([gt_array, gt_array2], axis=1)
    full_geno_matrix = pd.DataFrame(full_geno_matrix, dtype="float")
    corr_matrix = full_geno_matrix.corr()#.fillna(0)
    notnan_ids = list(corr_matrix.columns[np.logical_not(corr_matrix.isna().all().values)])
    corr_matrix = corr_matrix.loc[notnan_ids, notnan_ids]
    corr_matrix_array = corr_matrix.values
    np.savetxt(ldfile, corr_matrix_array, fmt='%.3f')
    print(corr_matrix.shape)

    PrintLine("Writing Z-score file", outf)
    pgc_snps_pval = snp_assoc_result
    pgc_snps = set(snps_rsid)

    pgc_strs_pval = str_assoc_result
    pgc_strs = set(str_rsid)

    pgc_strs_pval["sign"] = np.sign(np.log(pgc_strs_pval['OR']))
    pgc_strs_pval['Zscore'] = pgc_strs_pval["sign"]*abs(norm.ppf(pgc_strs_pval["P"]/2))
    pgc_strs_pval = pgc_strs_pval[pgc_strs_pval.ID.isin(pgc_strs)]
    pgc_strs_pval = pgc_strs_pval[['BP','ID','Zscore','P','vtype']]

    pgc_snps_pval = pgc_snps_pval[pgc_snps_pval.ID.isin(pgc_snps)]
    pgc_snps_pval = pgc_snps_pval[['BP','ID','Zscore','P','vtype']]

    str_snp = pd.concat([pgc_strs_pval,pgc_snps_pval])
    str_snp_caviar = str_snp.sort_values(by=['BP'])[['ID','Zscore']]
    str_snp_caviar = str_snp_caviar.set_index("ID")
    str_snp_caviar = str_snp_caviar.loc[list(corr_matrix.index)].reset_index()
    str_snp_caviar.to_csv(zfile, sep="\t", index=False, header=False)


    PrintLine("Running CAVIAR", outf)
    cmd = "/home/gymrek/ssaini/caviar/caviar/CAVIAR-C++/CAVIAR -o %s -l %s -z %s -c %s"%(outfile, ldfile, zfile, numcausal)
    p = Popen(cmd, shell=True, stdout=DEVNULL, stderr=DEVNULL)
    output = p.communicate()[0]
    if p.returncode != 0:
        PrintLine("CAVIAR failed", outf)
        return False

    caviar_output = pd.read_csv(outfile+"_post", delim_whitespace=True)
    caviar_output = caviar_output.merge(str_snp, left_on="SNP_ID", right_on="ID", how="left")
    caviar_output = caviar_output.sort_values(by=['Causal_Post._Prob.'], ascending=False)
    caviar_output = caviar_output.reset_index()
    caviar_output['index'] = caviar_output.index

    print("============ Top Causal Variants ============")
    print(caviar_output.head())

    print("============ Top STR Variants ============")
    print(caviar_output[caviar_output['vtype']=="STR"].head())

    fig_df = caviar_output.copy()
    fig_df = fig_df[['index','BP','P','vtype']]
    fig_df.loc[:,'logP'] = -np.log10(fig_df.P.values)

    sns.set(style="darkgrid")
    fig, ax = plt.subplots(figsize=(10,5))
    ax = sns.scatterplot(x="BP", y="logP", hue="vtype", data=fig_df, ax=ax, palette=dict(SNP="blue", STR="red"))
    ax.legend(bbox_to_anchor=(1, 1), ncol=1)
    xlabels = ['{:,.2f}'.format(x) + 'M' for x in ax.get_xticks()/1000000]
    ax.set_xticklabels(xlabels)
    ax.set_xlabel("BP")
    ax.set_ylabel("-log10 P")
    ax.set_title("Manhattan Plot %d:%d-%d"%(CHROM,START,END))
    fig.savefig("output_figs/fig1_chr%d-%d-%d.png"%(CHROM,START,END))

    fig, ax = plt.subplots(figsize=(10,5))
    ax = sns.scatterplot(x="index", y="logP", hue="vtype", data=fig_df, ax=ax, palette=dict(SNP="blue", STR="red"))
    ax.legend(bbox_to_anchor=(1, 1), ncol=1)
    ax.set_xlabel("CAVIAR Rank")
    ax.set_ylabel("-log10 P")
    ax.set_title("CAVIAR Rank vs. GWAS P-value %d:%d-%d"%(CHROM,START,END))
    fig.savefig("output_figs/fig2_chr%d-%d-%d.png"%(CHROM,START,END))

if __name__ == "__main__":
    main()
