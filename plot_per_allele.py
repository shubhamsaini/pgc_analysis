#!/usr/bin/env python
"""

Author: Shubham Saini
shubhamsaini@eng.ucsd.edu

"""
import matplotlib
matplotlib.use('Agg')


import argparse
import sys

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

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
    parser.add_argument("--assoc", help="PlinkSTR ouput", required=True, type=str)
    parser.add_argument("--min-maf", help="Minimum allele frequency. Default: 0.0", required=False, type=float, default=0.0)
    parser.add_argument("--out-dir", help="Output Directory", required=False, type=str, default=".")
    args = parser.parse_args()

    data = pd.read_csv(args.assoc, delim_whitespace=True)

    data = data.dropna()
    data = data[data['OR'] > 0]
    data = data[data['MAF'] > args.min_maf]

    data = data[data['SNP'].str.contains('length')]
    lengths = [int(x.split("-")[-1]) for x in data['SNP'].tolist()]
    pvals = [-np.log10(x) for x in data['P'].tolist()]
    odds = [x for x in data['OR'].tolist()]
    ci_low = [float(x.split("-")[0]) for x in data['CI95'].tolist()]
    ci_high = [float(x.split("-")[1]) for x in data['CI95'].tolist()]
    af = [round(x, 3) for x in data['MAF'].tolist()]
    plot_data = pd.DataFrame({'P':pvals, 'len':lengths, 'odds':odds, 'af':af, 'ci_low':ci_low, 'ci_high':ci_high})

    plot_data = plot_data.sort_values(by=['len'])

    fig = plt.figure()
    plt.errorbar(plot_data.len, plot_data.odds, yerr=(plot_data.odds-plot_data.ci_low), marker="o", linewidth=1, color="black")
    for line in range(0,plot_data.shape[0]):
        plt.text(plot_data.len[line], plot_data.odds[line]+((plot_data.odds[line]-plot_data.ci_low[line])/10), plot_data.af[line], horizontalalignment='left', color='black', weight='semibold')
    plt.title('chr%d:%d %s'%(data['CHROM'].values[0], data['BP'].values[0], data['SNP'].values[0].split("-")[0]))
    fig.savefig(args.out_dir+"/"+args.assoc+'.png')

if __name__ == "__main__":
    main()
