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
VCF=/mnt/scratch/wils1582/usa_capsella_noodles/cbp_all_msu_filt.vcf.gz
POPS=/mnt/home/wils1582/capsella_introgression/pixy_pops_NYC_ownpop.txt
OUTDIR=/mnt/scratch/wils1582
PREFIX=nyc_allSites_CBP

# Optional VARS
BED=/mnt/scratch/wils1582/degenotate_out/degeneracy-4-fold-sites.txt 

# # first, make sure the vcf is indexed
# tabix $VCF 
# 
 pixy --stats pi fst dxy \
	--sites_file $BED \
	--window_size 1 \
	--chromosomes 'jlSCF_10' \
 --vcf $VCF \
 --populations $POPS \
 --n_cores 10 \
 --output_folder "$OUTDIR" \
 --output_prefix "$PREFIX" \
	--bypass_invariant_check 'yes'
