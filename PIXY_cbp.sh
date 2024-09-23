#!/bin/bash
#
#SBATCH --job-name=PIXYcbp
#SBATCH --nodes=10
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-10:00:00
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wils1582@msu.edu
#SBATCH --output=/mnt/scratch/wils1582/slurm/slurm-%A.out
#
# Population statistics with PIXY on C. bursa-pastoris
# August 27, 2024

# PIXY population file generation in capsella_introgression_popgen.Rmd
cd /mnt/scratch/wils1582

# initialize conda environment
#conda activate pixy
source /mnt/home/wils1582/miniconda3/bin/activate pixy

#modules
module purge
module load tabixpp/1.1.2-GCC-12.3.0

#### VARS
VCF=/mnt/home/wils1582/allSites_CBP_final.filtered.vcf.gz
POPS=/mnt/home/wils1582/capsella_introgression/pixy_pops_NYC_ownpop.txt
OUTDIR=/mnt/scratch/wils1582
PREFIX=nyc_allSites_CBP

# Optional VARS
BED=/mnt/home/wils1582/capsella_introgression/nyc_rubella_pixy.bed

# # first, make sure the vcf is indexed
# tabix $VCF 
# 
# # Run PIXY on each chromosome in windows
# pixy --stats pi fst dxy \
# --vcf $VCF \
# --populations $POPS \
# --window_size 10000 \
# --n_cores 10 \
# --output_folder "$OUTDIR"/ \
# --output_prefix "$PREFIX"

# Run PIXY on each site (much slower than windows)
#pixy --stats pi fst dxy \
#--vcf $VCF \
#--populations $POPS \
#--window_size 100 \
#--n_cores 10 \
#--output_folder "$OUTDIR"/ \
#--output_prefix "$PREFIX" \
#--bypass_invariant_check 'yes' \
#--chromosomes 'jlSCF_10,jlSCF_11'

# # Run PIXY on intervals using a bed file
# pixy --stats pi fst dxy \
# --vcf $VCF \
# --populations $POPS \
# --bed_file "$BED" \
# --n_cores 10 \
# --output_folder "$OUTDIR" \
# --output_prefix "$PREFIX" \
# --bypass_invariant_check 'yes'
 pixy --stats pi fst dxy \
 --vcf $VCF \
 --populations $POPS \
 --window_size 10000 \
 --n_cores 10 \
 --output_folder "$OUTDIR" \
 --output_prefix "$PREFIX" \
 --bypass_invariant_check 'yes'
