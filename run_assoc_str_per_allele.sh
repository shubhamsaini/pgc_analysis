#!/bin/bash

###
### sbatch -t 0-12:00 --get-user-env ./run_assoc_str_locus_alleles.sh caviar_pcausal0.1.txt
###
###

cd /home/gymrek/ssaini/pipeline_per_locus/

LOCIFILE=$1

while read -r arg_1 arg_2 arg_3; do
    CHROM=${arg_1}
    POS=${arg_2}

    /home/gymrek/ssaini/gwas/plinkSTR_v4.py \
    --vcf /project/gymrek/chr${CHROM}/strs/pgc.strs.imputed.chr${CHROM}.vcf.gz \
    --covar /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/scz/wave2/v1/prune.bfile.cobg.PGC_SCZ49.sh2.menv.mds_cov \
    --covar-name C1,C2,C3,C4,C5,C6,C7,C9,C15,C18 \
    --out pgc.chr${CHROM}.${POS}.alleles.assoc  \
    --fam /home/gymrek/pgc_imputation/PGC_eur.fam \
    --sex --logistic --remove-rare-str-alleles 0.05 \
    --region ${CHROM}:${POS}-${POS} \
    --allele-tests --allele-tests-length

./plot_per_allele.py --assoc pgc.chr${CHROM}.${POS}.alleles.assoc
done < ${LOCIFILE}
