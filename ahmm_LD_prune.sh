#!/bin/bash
#
#SBATCH --job-name=ahmmLDPrune
#SBATCH --nodes=10
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-10:00:00
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wils1582@msu.edu
#SBATCH --output=/mnt/scratch/wils1582/slurm/slurm-%A_%a.out
#
# Ancestry HMM LD Pruning
# Prune LD SNPs in parental populations for Ancestry HMM analysis

# M. Wilson Brown
# July 3, 2024

### VARIABLES
VCF=/mnt/research/josephslab/Maya/capsella/vcf/adrian_vcf/final_filtered/NPCRCG_CBP_filtered_final.vcf.gz
OUTDIR=/mnt/home/wils1582/capsella_introgression
CR="$OUTDIR"/c_rubella.txt
AS_CBP="$OUTDIR"/eAsia_cbp.txt
ALL_SAMPLES="$OUTDIR"/sample_names.txt

# move directories
cd /mnt/scratch/wils1582/
# # purge modules
# module purge
# # Load modules
# ml PLINK/2.00a3.7-gfbf-2023a
# 
# # Sites in LD in C. rubella population
# plink2 --vcf $VCF \
# --indep-pairwise 100 5 0.2 \
# --keep $CR \
# --allow-extra-chr \
# --set-all-var-ids @:# \
# --out CR 
# 
# # Sites in LD in E. Asia C. bursa-pastoris population
# # NOTE: calculating LD on a small number of individuals is fraught, but alas
# plink2 --vcf $VCF \
# --indep-pairwise 100 5 0.2 \
# --keep $AS_CBP \
# --allow-extra-chr \
# --set-all-var-ids @:# \
# --bad-ld
# --out eAsia_CBP
# 
# # Combine sites
# cat eAsia_CBP.prune.in CR.prune.in > temp_keep_sites.prune.in

###### Missing variant frequency in parental pops
# Calculate missingness per variant in each parental population separately
plink2 --vcf "$VCF" \
  --keep $CR
  --extract temp_keep_sites.prune.in
  --missing variant-only vcols=chrom,pos,nmiss,nobs,fmiss \
  --allow-extra-chr \
  --double-id \
  --out cr_missing

# Same for parent2
plink2 --vcf "$VCF" \
  --keep $AS_CBP
  --extract temp_keep_sites.prune.in
  --missing variant-only vcols=chrom,pos,nmiss,nobs,fmiss \
  --allow-extra-chr \
  --double-id \
  --out as_missing

# Combine variant missing files; remove those with high missing frequency, write new sites to file
Rscript rmParentalMissing.R cr_missing.vmiss as_missing.vmiss 0.9

# Plink removes the allele depth feild from the vcf and I need that to convert to AHMM format
# so we will prune with bcftools (which I also forgot last time)

# remove them from kept sites file
# purge and load modules
module purge
ml BCFtools/1.19-GCC-13.2.0
# use bcftools to select only sites in keep_sites file
bcftools view --targets-file keep_filt_sites.txt $VCF \
  -Ov -o "$OUTDIR"/ahmm_pruned_all.vcf

# convert VCF to Ancestry HMM input format
python3 "$OUTDIR"/vcf2ahmm.py -v "$OUTDIR"/ahmm_pruned_all.vcf -s "$OUTDIR"/hmm_sample_mapping.txt
