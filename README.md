# pgc_analysis

# CAVIAR script to run on snorlax
```bash
./run_caviar_snorlax.py
--vcf /storage/s1saini/hipstr_allfilters/str_snp/tredparse/chr6.str.snp.with.tredparse.vcf.gz
--str-assoc /storage/s1saini/caviar/str_assoc/pgc.chr6.assoc
--snp-assoc /storage/s1saini/caviar/snp_assoc/chr6.snps.metal.txt
--snp-rsid /storage/s1saini/caviar/snp_assoc/chr6.pos.rsid.txt
--chrom 6 --pos 170870996
--window-kb 100
```

```bash
./run_caviar_lisa.py \
--str-vcf /nfs/gymrek/chr22/strs/pgc.strs.imputed.chr22.vcf.gz \
--snp-vcf /home/gymrek/ssaini/caviar/chr22.files.list \
--str-assoc /home/gymrek/ssaini/gwas/assoc_results/str_assoc/pgc.chr22.assoc \
--snp-assoc /home/gymrek/ssaini/gwas/assoc_results/snp_assoc/chr22.snps.metal.txt \
--snp-rsid /home/gymrek/ssaini/gwas/assoc_results/snp_assoc/chr22.pos.rsid.txt \
--chrom 22 --pos 39975317 \
--window-kb 10
```
