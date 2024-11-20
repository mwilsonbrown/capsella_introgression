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
VCF=/mnt/research/josephslab/Maya/capsella/vcf/filtered/CBP_CRCGCONP_maf_final.vcf.gz
OUTDIR=/mnt/home/wils1582/capsella_introgression
CR="$OUTDIR"/c_rubella.txt
AS_CBP="$OUTDIR"/eAsia_cbp.txt
ALL_SAMPLES="$OUTDIR"/sample_names.txt

# move directories
cd /mnt/scratch/wils1582/ahmm_workflow/
# # purge modules
module purge
# # Load modules
ml PLINK/2.00a3.7-gfbf-2023a
#
# # Sites in LD in C. rubella population
 plink2 --vcf $VCF \
 --indep-pairwise 500 50 0.5 \
 --keep $CR \
 --allow-extra-chr \
 --set-all-var-ids @:# \
 --out CR

 # Sites in LD in E. Asia C. bursa-pastoris population
 # NOTE: calculating LD on a small number of individuals is fraught, but alas
 plink2 --vcf $VCF \
 --indep-pairwise 500 50 0.5 \
 --keep $AS_CBP \
 --allow-extra-chr \
 --set-all-var-ids @:# \
 --bad-ld
 --out eAsia_CBP

# Combine sites
cat eAsia_CBP.prune.in CR.prune.in > temp_keep_sites.prune.in
# remove duplicate lines
awk '!seen[$0]++' temp_keep_sites.prune.in > keep_sites_nodups.prune.in

# for some reason, I need to duplicate the columns for the samples file for --missing
awk '$2=$1' $CR > cr_double.txt
awk '$2=$1' $AS_CBP > as_double.txt

##### Missing variant frequency in parental pops
#Calculate missingness per variant in each parental population separately
plink2 --vcf "$VCF" \
  --keep cr_double.txt \
  --extract keep_sites_nodups.prune.in \
  --missing variant-only vcols=chrom,pos,nmiss,nobs,fmiss \
  --allow-extra-chr \
  --double-id \
  --set-all-var-ids @:# \
  --out cr_missing

# Same for parent2
plink2 --vcf "$VCF" \
  --keep as_double.txt \
  --extract keep_sites_nodups.prune.in \
  --missing variant-only vcols=chrom,pos,nmiss,nobs,fmiss \
  --allow-extra-chr \
  --double-id \
  --set-all-var-ids @:# \
  --out as_missing

# modules 
module purge
module load R/4.4.1-gfbf-2023b

# Combine variant missing files; remove those with high missing frequency, write new sites to file
Rscript rmParentalMissing.R cr_missing.vmiss as_missing.vmiss 0.9

# Plink removes the allele depth feild from the vcf and I need that to convert to AHMM format
# so we will prune with bcftools (which I also forgot last time)

# remove them from kept sites file
# purge and load modules
module purge
ml BCFtools/1.19-GCC-13.2.0
# use bcftools to select only sites in keep_sites file; also remove NP here
bcftools view --targets-file keep_filt_sites.txt --samples ^ERR2990308.sam $VCF \
  -Ov -o ahmm_pruned_all.vcf

# convert VCF to Ancestry HMM input format
python3 "$OUTDIR"/edit_vcf2ahmm.py -v ahmm_pruned_all.vcf -s "$OUTDIR"/hmm_sample_mapping.txt
