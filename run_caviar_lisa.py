#!/usr/bin/env python3
"""

Author: Shubham Saini
shubhamsaini@eng.ucsd.edu


Run CAVIAR on SNP-STR GWAS results

Example usage:
./run_caviar_lisa.py \
--str-vcf /storage/s1saini/manuscript_strsnp/fig3/hipstr.1kg/eur/chr22/hipstr.chr22.eur.vcf.gz\
--snp-vcf /storage/s1saini/caviar/snp_gt_files.list \
--str-assoc /storage/s1saini/caviar/str_assoc/pgc.chr22.assoc \
--snp-assoc /storage/s1saini/caviar/snp_assoc/chr22.snps.metal.txt \
--snp-rsid /storage/s1saini/caviar/snp_assoc/chr22.pos.rsid.txt \
--chrom 22 --pos 39975317 \
--window-kb 100

"""


import argparse
import sys

from cyvcf2 import VCF
import numpy as np
import pandas as pd
from scipy.stats import norm
from subprocess import Popen, PIPE, DEVNULL
from collections import defaultdict

def PrintLine(text, f):
    f.write(text+"\n")
    f.flush()

def PROGRESS(msg, printit=True):
    if printit: # false for some messages when not in debug mode
        sys.stderr.write("%s\n"%msg.strip())

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
    ldfile = 'corr_mat_%d.txt'%(seed)
    zfile = 'zscore_%d.csv'%(seed)
    numcausal = "2"
    outfile = "caviar.output_%d"%(seed)

    PROGRESS('Seed %d'%(seed))
    PROGRESS("Loading STR GWAS results")

    str_assoc_result = pd.read_csv(args.str_assoc, delim_whitespace=True)
    str_assoc_result = str_assoc_result[[str_assoc_result.columns.tolist()[i] for i in [0,1,3,2,6,4]]]
    str_assoc_result.columns = ["CHR", "BP","P","ID","CI95","OR"]
    str_assoc_result = str_assoc_result[(str_assoc_result['BP']>=START) & (str_assoc_result['BP']<=END)]
    str_assoc_result = str_assoc_result[str_assoc_result.P != "NAN"]
    str_assoc_result['vtype'] = "STR"
    str_assoc_result['P'] = str_assoc_result.P.astype('float')
    str_assoc_result = str_assoc_result.sort_values(by=['BP'])

    PROGRESS("Loading SNP GWAS results")

    snp_assoc_result = pd.read_csv(args.snp_assoc, delim_whitespace=True)
    snp_pos_rsid = pd.read_csv(args.snp_rsid, delim_whitespace=True, names=["POS","MarkerName"])

    snp_assoc_result = snp_assoc_result.merge(snp_pos_rsid, how="inner", on="MarkerName")
    snp_assoc_result = snp_assoc_result[["POS","P-value","MarkerName","Zscore"]]
    snp_assoc_result.columns = ["BP","P","ID","Zscore"]
    snp_assoc_result['BP'] = snp_assoc_result.BP.astype('int')
    snp_assoc_result = snp_assoc_result[(snp_assoc_result['BP']>=START) & (snp_assoc_result['BP']<=END)]
    snp_assoc_result['vtype'] = "SNP"
    snp_assoc_result['CHR'] = CHROM
    snp_assoc_result['P'] = snp_assoc_result.P.astype('float')
    snp_assoc_result = snp_assoc_result.sort_values(by=['BP'])

    gwas_rsid = list(snp_assoc_result['ID'].values) + list(str_assoc_result['ID'].values)


    ### STR samples
    vcf = VCF(str_vcf_file, gts012=True)
    str_samples = vcf.samples

    PROGRESS("Loading SNP genotypes")
    gt_array = defaultdict(list)
    snps_rsid = set()
    str_rsid = set()
    samples = []

    for vcf_file_loc in snp_vcf_files:
        vcf = VCF(vcf_file_loc, gts012=True, samples = str_samples)
        samples = samples + list(vcf.samples)
        for v in vcf('%d:%d-%d'%(CHROM,START,END)):
            if len(v.REF) == 1 and str(v.ID) in gwas_rsid:
                gt = np.array(v.gt_types, dtype="float")
                gt[gt==3] = np.nan
                gt_array[str(v.ID)] = gt_array[str(v.ID)] + list(gt)
                snps_rsid.add(str(v.ID))

    gt_array = pd.DataFrame(gt_array, dtype="float")

    PROGRESS("Loading STR genotypes")
    gt_array2 = defaultdict(list)

    vcf = VCF(str_vcf_file, gts012=True, samples=samples)
    for v in vcf('%d:%d-%d'%(CHROM,START,END)):
        if str(v.ID) in gwas_rsid:
            gt = np.array(v.gt_types, dtype="float")
            gt[gt==3] = np.nan
            gt_array2[str(v.ID)] = gt_array2[str(v.ID)] + list(gt)
            str_rsid.add(str(v.ID))

    gt_array2 = pd.DataFrame(gt_array2, dtype="float")

    PROGRESS("Writing LD file")
    full_geno_matrix = pd.concat([gt_array, gt_array2], axis=1)
    full_geno_matrix = pd.DataFrame(full_geno_matrix, dtype="float")
    corr_matrix = full_geno_matrix.corr().fillna(0)**2
    corr_matrix_array = corr_matrix.values
    np.savetxt(ldfile, corr_matrix_array, fmt='%.3f')
    print(corr_matrix.shape)

    PROGRESS("Writing Z-score file")
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


    PROGRESS("Running CAVIAR")
    cmd = "/storage/s1saini/caviar/caviar/CAVIAR-C++/CAVIAR -o %s -l %s -z %s -c %s"%(outfile, ldfile, zfile, numcausal)
    p = Popen(cmd, shell=True, stdout=DEVNULL, stderr=DEVNULL)
    output = p.communicate()[0]
    if p.returncode != 0:
        PROGRESS("CAVIAR failed")
        return False

    caviar_output = pd.read_csv("caviar.output_%d_post"%(seed), delim_whitespace=True)
    caviar_output = caviar_output.merge(str_snp, left_on="SNP_ID", right_on="ID", how="left")
    caviar_output = caviar_output.sort_values(by=['Causal_Post._Prob.'], ascending=False)
    caviar_output = caviar_output.reset_index()
    caviar_output['index'] = caviar_output.index

    print(caviar_output.head())
    print(caviar_output[caviar_output['vtype']=="STR"].head())

if __name__ == "__main__":
    main()
