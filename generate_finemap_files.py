#!/usr/bin/env python
"""

Author: Shubham Saini
shubhamsaini@eng.ucsd.edu


Run FineMap on SNP-STR GWAS results

Example usage:
./generate_finemap_files.py \
--str-vcf /nfs/gymrek/chr22/strs/pgc.strs.imputed.chr22.vcf.gz \
--snp-vcf /home/gymrek/ssaini/caviar/chr22.files.list \
--str-assoc /home/gymrek/ssaini/gwas/str_assoc_feb2021/PGC.gp.iqr.filtered.assoc \
--snp-assoc /project/gymrek/snps_mega/assoc/PGC.maf_filt.assoc \
--snp-afreq /project/gymrek/snps_mega/afreq/PGC.maf_filt.afreq \
--chrom 22 --pos 39975317 \
--window-kb 10 \
--maxp 0.0001 \
--min-snps 100

"""


import argparse
import sys

from cyvcf2 import VCF
import numpy as np
import pandas as pd
from scipy.stats import norm
from collections import defaultdict
from itertools import combinations_with_replacement
from collections import Counter

from subprocess import Popen, PIPE, STDOUT

try:
    from subprocess import DEVNULL # py3k
except ImportError:
    import os
    DEVNULL = open(os.devnull, 'wb')

def PrintLine(text, f):
    f.write(text+"\n")
    f.flush()

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
    parser.add_argument("--snp-assoc", help="SNP GWAS file from Plink", required=True, type=str)
    parser.add_argument("--snp-afreq", help="SNP Allele Frequency File", required=True, type=str)
    parser.add_argument("--chrom", help="CHROM", required=True, type=str)
    parser.add_argument("--pos", help="POSITION", required=True, type=str)
    parser.add_argument("--window-kb", help="Window in KB", required=True, type=str)
    parser.add_argument("--maxp", help="Maximum p-value", required=False, type=float, default=0.0001)
    parser.add_argument("--min-snps", help="Minimum number of SNPs to include", required=False, type=int)
    parser.add_argument("--out-dir", help="Output directory name", required=False, type=str, default=".")
    parser.add_argument("--use-gp", help="Use GP field from Beagle output", action="store_true")
    parser.add_argument("--iqr-outliers", help="Filter outliers using IQR (GP-based regression only)", action="store_true")
    parser.add_argument("--iqr-outliers-min-samples", help="Minimum number of samples allowed per dosage. Use -1 to disable this.", default=100, type=int)

    parser.add_argument("--out", help="Output file name. Default: stdout", required=False, type=str)
    args = parser.parse_args()

    # Prepare output
    if args.out: outf = open(args.out,"w")
    else: outf = sys.stdout

    snp_vcf_files = []
    with open(args.snp_vcf) as f:
        for line in f:
            if len(line)>0:
                snp_vcf_files.append(line.strip())

    str_vcf_file = args.str_vcf

    CHROM = int(args.chrom)
    POS = int(args.pos)
    window=int(args.window_kb)*1000
    START = POS-window
    END = POS+window

    seed = np.random.randint(100000)
    ldfile = '%s/data.ld'%(args.out_dir)
    zfile = '%s/data.z'%(args.out_dir)
    outfile = "%s/caviar.output_%d_%d"%(args.out_dir,CHROM, POS)

    PrintLine('%d %d %d'%(CHROM, START, END), outf)
    PrintLine("Loading STR GWAS results", outf)

    str_assoc_result = pd.read_csv(args.str_assoc, delim_whitespace=True)
    str_assoc_result = str_assoc_result[["SNP", "CHROM", "BP", "ALLELE1", "ALLELE2", "MAF", "OR", "SE", "P"]]
    str_assoc_result = str_assoc_result[str_assoc_result['CHROM']==CHROM]
    str_assoc_result = str_assoc_result[(str_assoc_result['BP']>=START) & (str_assoc_result['BP']<=END)]
    str_assoc_result = str_assoc_result[str_assoc_result.P != "NAN"]
    str_assoc_result['vtype'] = "STR"
    str_assoc_result['P'] = str_assoc_result.P.astype('float')

    if args.maxp:
        str_assoc_result = str_assoc_result[str_assoc_result['P']<=float(args.maxp)]
    str_assoc_result = str_assoc_result.sort_values(by=['BP'])
    print(str_assoc_result.shape)
    if str_assoc_result.shape[0]==0:
        sys.exit('No STRs pass P-Val threshold')

    PrintLine("Loading SNP GWAS results", outf)

    snp_assoc_result = pd.read_csv(args.snp_assoc, delim_whitespace=True)
    snp_assoc_result = snp_assoc_result[snp_assoc_result['TEST']=="ADD"]
    snp_assoc_result = snp_assoc_result[snp_assoc_result['#CHROM']==CHROM]
    snp_assoc_result = snp_assoc_result[["ID", "#CHROM", "POS", "REF", "ALT", "OR", "SE", "P"]]
    snp_assoc_result.columns = ["SNP", "CHROM", "BP", "ALLELE1", "ALLELE2", "OR", "SE", "P"]
    snp_assoc_result['BP'] = snp_assoc_result.BP.astype('int')
    snp_assoc_result = snp_assoc_result[(snp_assoc_result['BP']>=START) & (snp_assoc_result['BP']<=END)]
    snp_assoc_result['vtype'] = "SNP"
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

    gwas_rsid = list(snp_assoc_result['SNP'].values) + list(str_assoc_result['SNP'].values)


    ### STR samples
    vcf = VCF(str_vcf_file)
    str_samples = vcf.samples


    region = '%d:%d-%d'%(CHROM,START,END)

    for vcf_file_loc in snp_vcf_files:
        cmd = f"bcftools query -r {region} -f '%ID\n' {vcf_file_loc} >> {args.out_dir}/snp_files_rsid_{CHROM}_{POS}.txt"
        print(cmd)
        p = Popen(cmd, shell=True, stdout=DEVNULL, stderr=DEVNULL)
        output = p.communicate()[0]
        if p.returncode != 0:
            print("Failed")

    snp_all_rsid = set()
    with open("%s/snp_files_rsid_%d_%d.txt"%(args.out_dir,CHROM, POS)) as f:
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
        vcf = VCF(vcf_file_loc, samples = str_samples)
        for v in vcf('%d:%d-%d'%(CHROM,START,END)):
            vcf_rsid.add(str(v.ID))
            if str(v.ID) in gwas_rsid and len(v.ALT)==1:
                gt = [i[0]+i[1] if (i[0]!=-1 and len(i)==3) else np.nan for i in v.genotypes]
                gt_array[str(v.ID)] = gt_array[str(v.ID)] + gt
                snps_rsid.add(str(v.ID))

        samples = samples + list(vcf.samples)
        not_found_rsid = snp_all_rsid.difference(vcf_rsid)
        for rsid in not_found_rsid:
            gt_array[rsid] = gt_array[rsid] + [np.nan]*len(vcf.samples)

    for i in list(gt_array):
        if len(gt_array[i]) == 0:
            del gt_array[i]

    gt_array = pd.DataFrame(gt_array, dtype="float")
    gt_array['sample'] = samples
    print(gt_array.shape)

    PrintLine("Loading STR genotypes",outf)
    gt_array2 = defaultdict(list)
    str_samples = []

    vcf = VCF(str_vcf_file, samples=samples)
    for v in vcf('%d:%d-%d'%(CHROM,START,END)):
        if str(v.ID) in gwas_rsid:
            if args.use_gp:
                genotypes = np.array(v.genotypes)[:,0:2]
                alleles_lengths = [len(v.REF)-len(v.REF)]+[len(i)-len(v.REF) for i in v.ALT]
                comb = list(combinations_with_replacement(range(len(alleles_lengths)), 2))
                idx = [gp_position(genotype[0], genotype[1]) for genotype in comb]
                geno_sum_lengths = [(alleles_lengths[genotype[0]] + alleles_lengths[genotype[1]]) for genotype in comb]
                #geno_sum_lengths_sorted = [geno_sum_lengths[i] for i in idx]
                geno_sum_lengths_dict = dict(zip(idx, geno_sum_lengths))
                geno_sum_lengths_sorted = [geno_sum_lengths_dict[i] for i in sorted(geno_sum_lengths_dict)]
                gp_sum = np.dot(geno_sum_lengths_sorted, v.format('GP').T)
                gt = list(gp_sum)

                if args.iqr_outliers:
                    q75, q25 = np.percentile(gp_sum, [75 ,25])
                    iqr = q75 - q25
                    low_threshold = q25 - (1.5*iqr)
                    high_threshold = q75 + (1.5*iqr)

                    if args.iqr_outliers_min_samples > 0:
                        alleles_idx = list(range(len(alleles_lengths)))
                        genotypes_lenghts = replace_allele_idx(genotypes, alleles_idx, alleles_lengths)
                        genotypes_lenghts = np.sum(genotypes_lenghts, axis=1)
                        gt_counts = Counter(genotypes_lenghts)
                        valid_gt_counts = [i for i in gt_counts if gt_counts[i] > args.iqr_outliers_min_samples]
                        if len(valid_gt_counts) > 0:
                            min_gt_counts = min(valid_gt_counts)
                            max_gt_counts = max(valid_gt_counts)
                            low_threshold = min(low_threshold, min_gt_counts)
                            high_threshold = max(high_threshold, max_gt_counts)

                    gt = [i if (i<=high_threshold and i>=low_threshold) else np.nan for i in gt]
            else:
                gt = [(len(i.split("|")[0]) + len(i.split("|")[1])) - 2*len(v.REF) if "." not in i else np.nan for i in v.gt_bases]
            gt_array2[str(v.ID)] = gt_array2[str(v.ID)] + gt
            str_rsid.add(str(v.ID))

    str_samples = str_samples + list(vcf.samples)
    gt_array2 = pd.DataFrame(gt_array2, dtype="float")
    gt_array2['sample'] = str_samples
    print(gt_array2.shape)

    PrintLine("Writing LD file", outf)
    #full_geno_matrix = pd.concat([gt_array, gt_array2], axis=1)
    full_geno_matrix = pd.merge(gt_array, gt_array2, on="sample")
    full_geno_matrix.to_csv("%s/geno_chr%d-%d-%d.csv"%(args.out_dir,CHROM,START,END))
    full_geno_matrix = full_geno_matrix.drop(['sample'], axis=1)
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

    pgc_strs_pval = pgc_strs_pval[pgc_strs_pval.SNP.isin(pgc_strs)]
    pgc_strs_pval['beta'] = np.log(pgc_strs_pval['OR'])
    pgc_strs_pval["sign"] = np.sign(np.log(pgc_strs_pval['OR']))
    pgc_strs_pval['Z'] = pgc_strs_pval["sign"]*abs(norm.ppf(pgc_strs_pval["P"]/2))
    pgc_strs_pval['beta_se'] = pgc_strs_pval['beta'] / pgc_strs_pval['Z']
    pgc_strs_pval = pgc_strs_pval[["SNP", "CHROM", "BP", "ALLELE1", "ALLELE2", "MAF", "beta", "beta_se"]]
    pgc_strs_pval.columns = ["rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se"]
    pgc_strs_pval['flip'] = 1

    pgc_snps_pval = pgc_snps_pval[pgc_snps_pval.SNP.isin(pgc_snps)]
    pgc_snps_pval['beta'] = np.log(pgc_snps_pval['OR'])
    pgc_snps_pval["sign"] = np.sign(np.log(pgc_snps_pval['OR']))
    pgc_snps_pval['Z'] = pgc_snps_pval["sign"]*abs(norm.ppf(pgc_snps_pval["P"]/2))
    pgc_snps_pval['beta_se'] = pgc_snps_pval['beta'] / pgc_snps_pval['Z']
    snp_afreq = pd.read_csv(args.snp_afreq, delim_whitespace=True)
    pgc_snps_pval = pgc_snps_pval.merge(snp_afreq, left_on="SNP", right_on="ID")
    pgc_snps_pval = pgc_snps_pval[["SNP", "CHROM", "BP", "ALLELE1", "ALLELE2", "ALT_FREQS", "beta", "beta_se"]]
    pgc_snps_pval.columns = ["rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se"]
    pgc_snps_pval['flip'] = 0

    str_snp = pd.concat([pgc_strs_pval,pgc_snps_pval])
    str_snp_caviar = str_snp.sort_values(by=['position'])
    str_snp_caviar = str_snp_caviar.set_index("rsid")
    str_snp_caviar = str_snp_caviar.loc[list(corr_matrix.index)].reset_index()
    str_snp_caviar = str_snp_caviar[["rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se", "flip"]]
    str_snp_caviar.to_csv(zfile, sep=" ", index=False)

if __name__ == "__main__":
    main()
