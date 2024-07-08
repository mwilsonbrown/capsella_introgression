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
# purge modules
module purge
# Load modules
ml PLINK/2.00a3.7-gfbf-2023a

# Sites in LD in C. rubella population
plink2 --vcf $VCF \
--indep-pairwise 100 5 0.2 \
--keep $CR \
--allow-extra-chr \
--set-all-var-ids @:# \
--out CR 

# Sites in LD in E. Asia C. bursa-pastoris population
# NOTE: calculating LD on a small number of individuals is fraught, but alas
plink2 --vcf $VCF \
--indep-pairwise 100 5 0.2 \
--keep $AS_CBP \
--allow-extra-chr \
--set-all-var-ids @:# \
--bad-ld
--out eAsia_CBP

# Combine sites
cat eAsia_CBP.prune.in CR.prune.in > keep_sites.prune.in

# Prune from admixed individuals
# plink2 --vcf $VCF \
# --extract keep_sites.prune.in \
# --keep $ALL_SAMPLES \
# --recode vcf bgz \
# --allow-extra-chr \
# --set-all-var-ids @:# \
# --out "$OUTDIR"/ahmm_pruned_cbp

# I am pretty sure plink removes the allele depth feild from the vcf and I need that
# so we will prune with bcftools (which I also forgot last time)

# reformat keep_sites file to two columns for bcftools; replaces colon with tab
sed -i 's/:/\t/g' keep_sites.prune.in 

# use bcftools to select only sites in keep_sites file
bcftools view --targets-file keep_sites.prune.in $VCF -Ov -o "$OUTDIR"/ahmm_pruned_all

# convert VCF to Ancestry HMM input format
python3 "$OUTDIR"/vcf2ahmm.py -v "$OUTDIR"/ahmm_pruned_all.vcf -s "$OUTDIR"/hmm_sample_mapping.txt
