#!/bin/bash


### usage: ./run_pipeline_per_locus.sh  6 170870996 167870996 173870996 1500



CHROM=$1
POS=$2
START=$3
END=$4
CAV_WINDOW=$5


ls /home/gymrek/ssaini/vcf/chr${CHROM}/*chr${CHROM}.vcf.gz > chr${CHROM}.files.list
bcftools merge -m id -l chr${CHROM}.files.list -Oz -o pgc.chr${CHROM}.${POS}.snps.vcf.gz -r ${CHROM}:${START}-${END}
bcftools index -t pgc.chr${CHROM}.${POS}.snps.vcf.gz

/home/gymrek/ssaini/plink2 \
--vcf pgc.chr${CHROM}.${POS}.snps.vcf.gz \
--covar prune.bfile.cobg.PGC_SCZ49.sh2.menv.mds_cov_vcf \
--covar-name C1,C2,C3,C4,C5,C6,C7,C9,C15,C18 \
--logistic \
--chr ${CHROM} \
--out pgc.chr${CHROM}.${POS}.snps.mega.plink \
--ci 0.95 \
--fam pgc.famfile.vcf.fam

cat pgc.chr${CHROM}.${POS}.snps.mega.plink.PHENO1.glm.logistic | awk '{if (NR>=2 && $7=="ADD") print;  else if (NR<=1) print;}' | sponge pgc.chr${CHROM}.${POS}.snps.mega.plink.PHENO1.glm.logistic
rm pgc.chr${CHROM}.${POS}.snps.mega.plink.PHENO1.glm.logistic.id

ls pgc.chr${CHROM}.${POS}.snps.vcf.gz  > chr${CHROM}.files.list

/home/gymrek/ssaini/caviar/run_caviar_lisa_v2.py \
--str-vcf /project/gymrek/chr${CHROM}/strs/pgc.strs.imputed.chr${CHROM}.vcf.gz \
--snp-vcf chr${CHROM}.files.list \
--str-assoc /home/gymrek/ssaini/caviar/str_assoc/pgc.chr${CHROM}.assoc \
--snp-assoc pgc.chr${CHROM}.${POS}.snps.mega.plink.PHENO1.glm.logistic \
--snp-rsid /home/gymrek/ssaini/caviar/snp_assoc/chr${CHROM}.pos.rsid.txt \
--chrom ${CHROM} --pos ${POS} \
--window-kb ${CAV_WINDOW} \
--maxp 0.0001
