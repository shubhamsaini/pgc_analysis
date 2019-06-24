#!/usr/bin/env python3
"""

Author: Shubham Saini
shubhamsaini@eng.ucsd.edu


Run CAVIAR on SNP-STR GWAS results

Example usage:
./run_caviar_snorlax.py \
--vcf /storage/s1saini/hipstr_allfilters/str_snp/tredparse/chr6.str.snp.with.tredparse.vcf.gz \
--str-assoc /storage/s1saini/caviar/str_assoc/pgc.chr6.assoc \
--snp-assoc /storage/s1saini/caviar/snp_assoc/chr6.snps.metal.txt \
--snp-rsid /storage/s1saini/caviar/snp_assoc/chr6.pos.rsid.txt \
--chrom 6 --pos 170870996 \
--window-kb 100

"""


import argparse
import sys

from cyvcf2 import VCF
import numpy as np
import pandas as pd
from scipy.stats import norm
from subprocess import Popen, PIPE, DEVNULL

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
    parser.add_argument("--vcf", help="VCF File for LD calculation", required=True, type=str)
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
    
    input_vcf = args.vcf
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


    PROGRESS("Loading genotypes")
    vcf = VCF(input_vcf)
    gt_array = []
    snps_rsid = []
    str_rsid = []
    for v in vcf('%d:%d-%d'%(CHROM,START,END)):
        if len(v.REF) == 1:
            if str(v.ID) in gwas_rsid:
                gt_array.append(list(np.sum(np.array(v.genotypes)[:,0:2], axis=1)))
                snps_rsid.append(str(v.ID))
        else:
            str_gt = []
            if str(v.ID) in gwas_rsid:
                for gt in v.gt_bases:
                    str_gt.append(np.sum([len(i)-len(v.REF) for i in gt.split("|")]))
                gt_array.append(str_gt)
                str_rsid.append(str(v.ID))

    pgc_snps_pval = snp_assoc_result
    pgc_snps = set(snps_rsid)

    pgc_strs_pval = str_assoc_result
    pgc_strs = set(str_rsid)

    PROGRESS("Writing LD file")

    df = pd.DataFrame(np.transpose(gt_array))
    corr_matrix = df.corr().fillna(0)**2
    corr_matrix = corr_matrix.values
    np.savetxt(ldfile, corr_matrix, fmt='%.3f')
    corr_matrix.shape

    PROGRESS("Writing Z-score file")

    pgc_strs_pval["sign"] = np.sign(np.log(pgc_strs_pval['OR']))
    pgc_strs_pval['Zscore'] = pgc_strs_pval["sign"]*abs(norm.ppf(pgc_strs_pval["P"]/2)) # 6.1
    pgc_strs_pval = pgc_strs_pval[pgc_strs_pval.ID.isin(pgc_strs)]
    pgc_strs_pval = pgc_strs_pval[['BP','ID','Zscore','P','vtype']]

    pgc_snps_pval = pgc_snps_pval[pgc_snps_pval.ID.isin(pgc_snps)]
    pgc_snps_pval = pgc_snps_pval[['BP','ID','Zscore','P','vtype']]

    str_snp = pd.concat([pgc_strs_pval,pgc_snps_pval])
    str_snp_caviar = str_snp.sort_values(by=['BP'])[['ID','Zscore']]
    str_snp_caviar.to_csv(zfile, sep="\t", index=False, header=False)


    PROGRESS("Running CAVIAR")
    cmd = "/storage/s1saini/caviar/caviar/CAVIAR-C++/CAVIAR -o %s -l %s -z %s -c %s"%(outfile, ldfile, zfile, numcausal)
    p = Popen(cmd, shell=True, stdout=DEVNULL, stderr=DEVNULL)
    output = p.communicate()[0]
    if p.returncode != 0:
        PROGRESS("CAVIAR failed")
        return False

    caviar_output = pd.read_csv("caviar.output_post", delim_whitespace=True)
    caviar_output = caviar_output.merge(str_snp, left_on="SNP_ID", right_on="ID", how="left")
    caviar_output = caviar_output.sort_values(by=['Causal_Post._Prob.'], ascending=False)
    caviar_output = caviar_output.reset_index()
    caviar_output['index'] = caviar_output.index

    print(caviar_output.head())
    print(caviar_output[caviar_output['vtype']=="STR"].head())
        
if __name__ == "__main__":
    main()
