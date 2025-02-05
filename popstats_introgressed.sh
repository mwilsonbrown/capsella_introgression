#!/bin/bash --login
#
#SBATCH --job-name=Pi_introgressed
#SBATCH --nodes=10
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-15:00:00
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
#SBATCH --export=NONE
#SBATCH --mem-per-cpu=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wils1582@msu.edu
#SBATCH --array=1-16%3
#SBATCH --output=/mnt/scratch/wils1582/slurm/slurm-%A.out
#
# Population statistics on introgressed regions with PIXY on C. bursa-pastoris
# Feb 5, 2025

# PIXY population file generation in capsella_introgression_popgen.Rmd
cd /mnt/scratch/wils1582/popstats_intro

# initialize conda environment
#conda activate pixy
source /mnt/home/wils1582/miniconda3/bin/activate pixy

#### VARS
VCF=/mnt/research/josephslab/Maya/capsella/vcf/filtered/CBP_allSites_msu_maf_reheader.vcf.gz
POPS=/mnt/home/wils1582/capsella_introgression/pixy_pops_NYC_ownpop.txt
OUTDIR=/mnt/scratch/wils1582/popstats_intro
PREFIX=all_multiinter_popstats
# Set BED file as first 3 columns from multiinter output file (ommitting the first line which is a header)
BED=$(tail -n +2 /mnt/scratch/wils1582/multiinter_bed/all_multiinter_rubella.txt | cut -f1,2,3)

# generate output space
#mkdir -p $OUTDIR

# calculate pi within introgressed regions
 pixy --stats pi \
	--bed_file $BED \
	--vcf $VCF \
 	--populations $POPS \
 	--n_cores 10 \
 	--output_folder "$OUTDIR" \
 	--output_prefix "$PREFIX" \
