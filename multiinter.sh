#!/bin/bash
#
#SBATCH --job-name=MULTIinter
#SBATCH --nodes=10
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-8:00:00
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wils1582@msu.edu
#SBATCH --output=/mnt/scratch/wils1582/slurm/slurm-%A_%a.out
#
# BEDtools Multi-inter
# July 12, 2024
# Maya Wilson Brown
#
#
# move to correct directory
cd /mnt/scratch/wils1582/bed_files/

#####VARS
OUTDIR=/mnt/home/capsella_introgression

# load modules
module purge
module load BEDTools/2.31.0-GCC-12.3.0

# output list of overlap of segments
bedtools multiinter -header -i *NY_rubella.bed > "$OUTDIR"/ny_multiinter_rubella.txt
bedtools multiinter -header -i *N_Europe_rubella.bed > "$OUTDIR"/neu_multiinter_rubella.txt
bedtools multiinter -header -i *MENA_rubella.bed > "$OUTDIR"/mena_multiinter_rubella.txt
