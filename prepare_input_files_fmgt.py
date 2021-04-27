#!/usr/bin/env python
"""

Author: Shubham Saini
shubhamsaini@eng.ucsd.edu

"""

import argparse
import sys

from cyvcf2 import VCF
import numpy as np
import pandas as pd
from collections import defaultdict
from itertools import combinations_with_replacement

def gp_position(j, k):
    p1 = (k*(k+1)/2)+j
    p2 = (j*(j+1)/2)+k
    return int(max(p1, p2))

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
    parser.add_argument("--str-snp-vcf", help="SNP VCF files list for LD calculation", required=True, type=str)
    parser.add_argument("--str-assoc", help="STR GWAS file from plinkSTR", required=True, type=str)
    parser.add_argument("--snp-assoc", help="SNP GWAS file from Metal", required=True, type=str)
    parser.add_argument("--fam", help="FAM File", required=True, type=str)
    parser.add_argument("--cov", help="Covariates file", type=str)
    parser.add_argument("--covar-name", help="Names of covariates to load. Comma-separated", type=str)
    parser.add_argument("--vcf-samples-delim", help="FID and IID delimiter in VCF", type=str)
    parser.add_argument("--chrom", help="CHROM", required=True, type=str)
    parser.add_argument("--pos", help="POSITION", required=True, type=str)
    parser.add_argument("--window-kb", help="Window in KB", required=True, type=str)
    parser.add_argument("--maxp", help="Maximum p-value", required=False, type=str)
    parser.add_argument("--maxp-strs", help="Maximum p-value", required=False, type=str)
    parser.add_argument("--min-snps", help="Minimum number of SNPs to include", required=False, type=int)
    parser.add_argument("--use-gp", help="Use GP field from Beagle output", action="store_true")
    parser.add_argument("--out", help="Output file name. Default: stdout", required=False, type=str)
    args = parser.parse_args()

    if args.out: outf = open(args.out,"w")
    else: outf = sys.stdout

    snp_vcf_files = []
    with open(args.str_snp_vcf) as f:
        for line in f:
            snp_vcf_files.append(line.strip())

    CHROM = int(args.chrom)
    POS = int(args.pos)
    window=int(args.window_kb)*1000
    START = POS-window
    END = POS+window


    PrintLine('%d %d %d'%(CHROM, START, END), outf)
    PrintLine("Loading STR GWAS results", outf)

    str_assoc_result = pd.read_csv(args.str_assoc, delim_whitespace=True)
    str_assoc_result = str_assoc_result[["CHROM", "BP","SNP","P","CI95","OR"]]
    str_assoc_result = str_assoc_result[str_assoc_result['CHROM']==CHROM]
    str_assoc_result = str_assoc_result[(str_assoc_result['BP']>=START) & (str_assoc_result['BP']<=END)]
    str_assoc_result = str_assoc_result[str_assoc_result.P != "NAN"]
    str_assoc_result['vtype'] = "STR"
    str_assoc_result['P'] = str_assoc_result.P.astype('float')

    if args.maxp_strs:
        str_assoc_result = str_assoc_result[str_assoc_result['P']<=float(args.maxp_strs)]
    str_assoc_result = str_assoc_result.sort_values(by=['BP'])
    print(str_assoc_result.shape)
    if str_assoc_result.shape[0]==0:
        sys.exit('No STRs pass P-Val threshold')

    PrintLine("Loading SNP GWAS results", outf)

    snp_assoc_result = pd.read_csv(args.snp_assoc, delim_whitespace=True)
    snp_assoc_result = snp_assoc_result[snp_assoc_result['TEST']=="ADD"]
    snp_assoc_result = snp_assoc_result[snp_assoc_result['#CHROM']==CHROM]
    snp_assoc_result = snp_assoc_result[["POS","P","ID","Z_STAT"]]
    snp_assoc_result.columns = ["BP","P","ID","Zscore"]
    snp_assoc_result['BP'] = snp_assoc_result.BP.astype('int')
    snp_assoc_result = snp_assoc_result[(snp_assoc_result['BP']>=START) & (snp_assoc_result['BP']<=END)]
    snp_assoc_result['vtype'] = "SNP"
    snp_assoc_result['CHROM'] = CHROM
    snp_assoc_result['P'] = snp_assoc_result.P.astype('float')
    if args.maxp:
        if args.min_snps and snp_assoc_result[snp_assoc_result['P']<=float(args.maxp)].shape[0] < args.min_snps:
            snp_assoc_result = snp_assoc_result.sort_values(by=['P'])
            snp_assoc_result = snp_assoc_result.head(args.min_snps)
        else:
            snp_assoc_result = snp_assoc_result[snp_assoc_result['P']<=float(args.maxp)]
    snp_assoc_result = snp_assoc_result.sort_values(by=['BP'])
    print(snp_assoc_result.shape)
    if snp_assoc_result.shape[0]==0:
        sys.exit('No SNPs pass P-Val threshold')

    gwas_rsid = list(snp_assoc_result['ID'].values) + list(str_assoc_result['SNP'].values)

    region = '%d:%d-%d'%(CHROM,START,END)

    PrintLine("Loading SNP genotypes", outf)

    snps_rsid = set()
    samples = []

    gt_array = defaultdict(list)
    for rsid in gwas_rsid:
        gt_array[rsid] = []

    for vcf_file_loc in snp_vcf_files:
        print("Loading genotypes from ",vcf_file_loc)
        vcf_rsid = set()
        vcf = VCF(vcf_file_loc)
        for v in vcf('%d:%d-%d'%(CHROM,START,END)):
            if str(v.ID) in gwas_rsid:
                str_gt = []
                vcf_rsid.add(str(v.ID))
                if len(v.REF)==1 and len(v.ALT)==1 and len(v.ALT[0])==1:
                    str_gt = list(np.sum(np.array(v.genotypes)[:,0:2], axis=1))
                else:
                    if args.use_gp:
                        genotypes = np.array(v.genotypes)[:,0:2]
                        alleles_lengths = [len(v.REF)]+[len(i) for i in v.ALT]
                        comb = list(combinations_with_replacement(range(len(alleles_lengths)), 2))
                        idx = [gp_position(genotype[0], genotype[1]) for genotype in comb]
                        geno_sum_lengths = [(alleles_lengths[genotype[0]] + alleles_lengths[genotype[1]]) for genotype in comb]
                        #geno_sum_lengths_sorted = [geno_sum_lengths[i] for i in idx]
                        geno_sum_lengths_dict = dict(zip(idx, geno_sum_lengths))
                        geno_sum_lengths_sorted = [geno_sum_lengths_dict[i] for i in sorted(geno_sum_lengths_dict)]
                        gp_sum = np.dot(geno_sum_lengths_sorted, v.format('GP').T)
                        str_gt = list(gp_sum)
                    else:
                        str_gt = [(len(i.split("|")[0]) + len(i.split("|")[1])) - 2*len(v.REF) for i in v.gt_bases]
                gt_array[str(v.ID)] = gt_array[str(v.ID)] + list(str_gt)
                snps_rsid.add(str(v.ID))

        samples = samples + list(vcf.samples)
        not_found_rsid = set(gwas_rsid).difference(vcf_rsid)
        for rsid in not_found_rsid:
            gt_array[rsid] = gt_array[rsid] + [np.nan]*len(vcf.samples)

    for i in list(gt_array):
        if len(gt_array[i]) == 0:
            del gt_array[i]

    gt_array = pd.DataFrame(gt_array, dtype="float")
    print(gt_array.shape)

    PrintLine("Writing Genotypes file", outf)
    full_geno_matrix = gt_array
    full_geno_matrix = pd.DataFrame(full_geno_matrix, dtype="float")
    full_geno_matrix['sample'] = samples


    PrintLine("Loading phenotype file", outf)
    pheno = pd.read_csv(args.fam, delim_whitespace=True, names=["FID", "IID", "Father_ID", "Mother_ID", "sex", "phenotype"])
    pheno = pheno[pheno["phenotype"].apply(str) != -9]
    pheno["phenotype"] = pheno["phenotype"].apply(int)-1 # convert to 0/1
    if args.vcf_samples_delim is not None:
        pheno["sample"] = pheno.apply(lambda x: x["FID"]+args.vcf_samples_delim+x["IID"], 1)
    else:
        pheno["sample"] = pheno.apply(lambda x: x["IID"], 1)

    PrintLine("Loading covariate file", outf)
    default_cols = ["FID", "IID"]
    if args.covar_name:
        colnames = default_cols+args.covar_name.split(",")
        covarcols = args.covar_name.split(",")
        cov = pd.read_csv(args.cov, delim_whitespace=True, usecols=colnames)
    pheno = pd.merge(pheno, cov, on=["FID","IID"])

    rmcovars = []
    for col in covarcols:
        if np.var(pheno[col])==0: rmcovars.append(col)

    geno_samples = set(full_geno_matrix['sample'])
    full_geno_matrix = full_geno_matrix.set_index('sample')
    pheno_samples = set(pheno["sample"])
    pheno = pheno.set_index('sample')

    non_missing_samples = list(geno_samples.intersection(pheno_samples))
    print("Number of samples: ", len(non_missing_samples))

    full_geno_matrix = full_geno_matrix.loc[non_missing_samples]
    full_geno_matrix.dropna(axis=1, how="all", inplace=True)
    full_geno_matrix.to_csv("genotypes.%d_%d.txt"%(CHROM, POS), sep="\t", index=False, na_rep="NA")

    PrintLine("Writing phenotype file", outf)
    pheno = pheno.loc[non_missing_samples]
    pheno.reset_index(inplace=True)
    pheno.to_csv("pheno.%d_%d.txt"%(CHROM, POS), columns=['phenotype'], header=False, index=False)

    PrintLine("Writing covariates file", outf)
    pheno.to_csv("cov.%d_%d.txt"%(CHROM, POS), columns=[item for item in covarcols if item not in rmcovars], header=False, index=False)



if __name__ == "__main__":
    main()
