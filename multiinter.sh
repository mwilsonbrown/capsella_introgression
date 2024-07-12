#!/bin/bash








# move to correct directory
cd mnt/scratch/wils1582/bed_files/

# load modules
module purge
module load GCC/10.2.0
module load BEDTools/2.30.0

# output list of overlap of segments
bedtools multiinter -header -i *_rubella_chunks.bed > nyc_multiinter_rubella.txt