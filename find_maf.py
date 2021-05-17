#!/usr/bin/env python
"""

Author: Shubham Saini
shubhamsaini@eng.ucsd.edu


Usage: ./find_af.py \
--vcf /storage/resources/datasets/SSC_SNP_v3/shapeit.chr22.with.ref.v3.vcf.gz \
--region 22:16788134-17788134


"""
import numpy as np
import pandas as pd
from cyvcf2 import VCF
import allel

import argparse
import sys

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def calcMAF(counts):
    return np.sum(np.sort(counts)[:-1])


def isSTR(record_alleles_to_bases):
    if record_alleles_to_bases[0] == 1 and record_alleles_to_bases[1]==1 and len(record_alleles_to_bases.keys())==2:
        return False
    else:
        return True


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--vcf", help="Input VCF File", required=True, type=str)
    parser.add_argument("--region", help="Fetch only chr:start-end from the VCF", required=False, type=str)
    args = parser.parse_args()



    VCF_FILE = args.vcf

    if args.region:
        region = args.region
    else:
        region=None


### Load meta data Information per record
### Used for finding STR allele lengths
### Used to create bim_data DF when BIM file not used
    meta_data = []
    vcf = VCF(VCF_FILE)
    if args.region:
        vcf_iter = vcf(args.region)
    else:
        vcf_iter = vcf()
    alleles_to_bases = {}
    for v in vcf_iter:
        alleles = [len(v.REF)]+[len(i) for i in v.ALT]
        alleles_to_bases[str(v.ID)] = dict(zip(range(len(alleles)), alleles))
        meta_data.append([v.CHROM, v.ID, v.POS])

    bim_data = pd.DataFrame(meta_data, columns=["CHR","SNP","BP"])
    bim_snps = list(bim_data['SNP'].values)


### Load Genotypes from the VCF ###
    callset = allel.read_vcf(VCF_FILE, region=region)
    genotype_array = callset['calldata/GT']
    print("ID, MAF, ALLELE:AF")
    for i in range(genotype_array.shape[0]):
        ID = callset['variants/ID'][i]
        record_alleles_to_bases = alleles_to_bases[ID]
        is_str = isSTR(record_alleles_to_bases)
        if not is_str:
            numSamples = callset['samples'].shape[0]
            alleles, counts = np.unique(genotype_array[i,:,:], return_counts=True)
            counts = counts/(2*numSamples)
            maf = calcMAF(counts)
        else:
            numSamples = callset['samples'].shape[0]
            u,inv = np.unique(genotype_array[i,:,:],return_inverse = True) # map allele index to allele length
            filtered_gt = np.array([record_alleles_to_bases.get(x,x) for x in u])[inv].reshape(genotype_array[i,:,:].shape) # map allele index to allele length
            alleles, counts = np.unique(filtered_gt, return_counts=True)
            counts = counts/(2*numSamples)
            maf = calcMAF(counts)
    output = " ".join(["%d:%f"%(i,j) for i,j in zip(alleles,counts)])
    print(ID, str(maf), output)




if __name__ == "__main__":
    main()
