#!/usr/bin/env python
"""

Author: Shubham Saini
shubhamsaini@eng.ucsd.edu

Usage: ./filter_gp.py \
--input-vcf /storage/s1saini/hipstr_allfilters/phased_feb18/hipstr.chr22.phased.vcf.gz

"""


import argparse
import sys

from cyvcf2 import VCF
import re
import numpy as np

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def write_variant(writer, chrom, pos, rsid, ref, alt, qual, filt, info, fmt, gt):
    fmt_variant = [str(chrom), str(pos), str(rsid), str(ref), str(alt), qual, filt, info, fmt] + gt
    writer.write("\t".join(fmt_variant)+"\n")

def gp_position(j,k):
    p1 = (k*(k+1)/2)+j
    p2 = (j*(j+1)/2)+k
    return int(max(p1,p2))

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--input-vcf", help="Input VCF", required=True, type=str)
    parser.add_argument("--output-vcf", help="Output VCF", required=False, type=str)
    parser.add_argument("--region", help="Region (optional)", required=False, type=str)
    parser.add_argument("--min-gp", help="Minimum GP", required=False, type=float, default=0.8)
    args = parser.parse_args()

    if args.output_vcf:
        output_file = args.output_vcf
        writer = open(output_file, "w")
    else:
        writer = sys.stdout
    writer.write("##fileformat=VCFv4.1\n")
    writer.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    writer.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    for i in range(1,23):
        writer.write('##contig=<ID=%d>\n'%i)

    VCF_FILE = args.input_vcf
    vcf = VCF(VCF_FILE)
    samples = vcf.samples

    header_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + samples
    writer.write("\t".join(header_fields)+"\n")

    vcf = VCF(VCF_FILE)
    if args.region:
        region = args.region
        for v in vcf(region):
            gtlist = v.genotypes
            gt_prob = v.format("GP")
            if len(gtlist) != len(gt_prob):
                continue
            str_gt = []
            for sample in range(len(gtlist)):
                gt_sample = gtlist[sample][0:2]
                gt_prob_sample = gt_prob[sample][gp_position(gt_sample[0], gt_sample[1])]
                if gt_prob_sample < args.min_gp:
                    gtlist[sample][0], gtlist[sample][1] = ".","."
                    str_gt.append("|".join([gtlist[sample][0], gtlist[sample][1]]))
                else:
                    str_gt.append("%d|%d"%(gtlist[sample][0], gtlist[sample][1]))
            write_variant(writer, v.CHROM, v.POS, v.ID, v.REF, ",".join(v.ALT), ".", "PASS", ".", "GT", str_gt)
    else:
        for v in vcf():
            gtlist = v.genotypes
            gt_prob = v.format("GP")
            if len(gtlist) != len(gt_prob):
                continue
            str_gt = []
            for sample in range(len(gtlist)):
                gt_sample = gtlist[sample][0:2]
                gt_prob_sample = gt_prob[sample][gp_position(gt_sample[0], gt_sample[1])]
                if gt_prob_sample < args.min_gp:
                    gtlist[sample][0], gtlist[sample][1] = ".","."
                    str_gt.append("|".join([gtlist[sample][0], gtlist[sample][1]]))
                else:
                    str_gt.append("%d|%d"%(gtlist[sample][0], gtlist[sample][1]))
            write_variant(writer, v.CHROM, v.POS, v.ID, v.REF, ",".join(v.ALT), ".", "PASS", ".", "GT", str_gt)

    writer.close()
if __name__ == "__main__":
    main()
