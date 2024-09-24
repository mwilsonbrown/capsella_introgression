# Test Bed file extraction from VCF for use in LOSTRUCT
# M. Wilson Brown
# September 23, 2024

# This extracts regions inferred to have C.rubella introgression from a VCF

# working directory
cd /mnt/scratch/wils1582/

# load bcftools
module purge
module load BCFtools/1.19-GCC-13.2.0

# NYC Capsella bursa-pastoris C. rubella regions
cut -f 1,2,3 /mnt/home/wils1582/capsella_introgression/nyc_multiinter_rubella.txt > ./nyc_cbp.bed
sed -i '1d' nyc_cbp.bed

# BCFtools to extract SNPs in introgressed regions
# sed command recodes genotypes as numeric
# because regions file had '.bed' suffix, regions are considered 0-based, half open
# exclude Neslia paniculata, but leave in all other samples.

bcftools query -f '%CHROM\t%POS[\t%GT]\n' /mnt/home/wils1582/ahmm_pruned_all.vcf.gz \
  --regions-file nyc_cbp.bed \
  --samples ^'ERR2990308.sam' |sed -e 's_|_/_g' > nyc_rubella_unphased_vcf.txt

# use sed to convert to numeric genotype calls  
sed -e 's_0/0_0_g' -e 's_\(0/1\)\|\(1/0\)_1_g' -e 's_1/1_2_g' -e 's_./._NA_g' nyc_rubella_unphased_vcf.txt > nyc_rubella_vcf.txt
