#!/bin/bash
#
#SBATCH --job-name=Pi4_bed_1based
#SBATCH --nodes=10
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-10:00:00
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
#SBATCH --export=NONE
#SBATCH --mem-per-cpu=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wils1582@msu.edu
#SBATCH --array=1-16%3
#SBATCH --output=/mnt/scratch/wils1582/slurm/slurm-%A_%a.out
#
# Population statistics with PIXY on C. bursa-pastoris
# August 27, 2024

# PIXY population file generation in capsella_introgression_popgen.Rmd
cd /mnt/scratch/wils1582

# initialize conda environment
#conda activate pixy
source /mnt/home/wils1582/miniconda3/bin/activate pixy

#### VARS
VCF=/mnt/research/josephslab/Maya/capsella/vcf/filtered/CBP_allsites_msu_maf_SNPsOnly.vcf.gz
POPS=/mnt/home/wils1582/capsella_introgression/pixy_pops_NYC_ownpop.txt
OUTDIR=/mnt/scratch/wils1582/4-fold_pi_1-based
PREFIX=nyc_allSites_CBP
CHROM=$(sed -n "$SLURM_ARRAY_TASK_ID"p /mnt/scratch/wils1582/degenotate_out/scaffold_names.txt)

# Optional VARS
#SITES=/mnt/scratch/wils1582/degenotate_out/degeneracy-4-fold-sites.txt 
#
# PIXY expects 1-indexed, not 0-indexed bed files so I created that file like so
# awk -d\t '$2=$2+1, $3=$3+1' degeneracy-4-fold_test.bed > degeneracy-4-fold-sites-1-based.bed
# where test.bed is just a copy of degeneracy-4-fold_sites.bed
BED=/mnt/scratch/wils1582/degenotate_out/degeneracy-4-fold-sites-1-based.bed

# generate output space
mkdir -p $OUTDIR

# # calculate pi at 4-fold degenerate sites
#  pixy --stats pi \
# 	--sites_file $SITES \
# 	--window_size 1 \
# 	--chromosomes "$CHROM" \
#  --vcf $VCF \
#  --populations $POPS \
#  --n_cores 10 \
#  --output_folder "$OUTDIR" \
#  --output_prefix "$PREFIX"_"$CHROM" \
# 	--bypass_invariant_check 'yes'

# calculate pi at 4-fold degenerate sites
 pixy --stats pi \
	--bed_file $BED \
	--chromosomes "$CHROM" \
 --vcf $VCF \
 --populations $POPS \
 --n_cores 10 \
 --output_folder "$OUTDIR" \
 --output_prefix "$PREFIX"_"$CHROM" \
	--bypass_invariant_check 'yes'
