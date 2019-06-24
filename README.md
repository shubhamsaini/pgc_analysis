# pgc_analysis

# CAVIAR script to run on snorlax
  ./run_caviar_snorlax.py \
  --vcf /storage/s1saini/hipstr_allfilters/str_snp/tredparse/chr6.str.snp.with.tredparse.vcf.gz \
  --str-assoc /storage/s1saini/caviar/str_assoc/pgc.chr6.assoc \
  --snp-assoc /storage/s1saini/caviar/snp_assoc/chr6.snps.metal.txt \
  --snp-rsid /storage/s1saini/caviar/snp_assoc/chr6.pos.rsid.txt \
  --chrom 6 --pos 170870996 \
  --window-kb 100
