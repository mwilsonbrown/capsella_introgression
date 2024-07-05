
# Ancestry HMM LD Pruning
# Prune LD SNPs in parental populations for Ancestry HMM analysis

# M. Wilson Brown
# July 3, 2024


### VARIABLES
VCF=
OUTDIR=
CR=
AS_CBP=
ALL_SAMPLES=

# purge modules
module purge
# Load modules
ml PLINK/2.00a3.7-gfbf-2023a

# 

# Sites in LD in C. rubella population
plink2 --vcf $VCF \
--indep-pairwise 100 5 0.2 \
--keep $CR \
--allow-extra-chr \
--set-all-var-ids @:# \
--out CR 

# Sites in LD in E. Asia C. bursa-pastoris population
plink2 --vcf $VCF \
--indep-pairwise 100 5 0.2 \
--keep $AS_CBP \
--allow-extra-chr \
--set-all-var-ids @:# \
--out eAsia_CBP

# Combine sites
cat eAsia_CBP.prune.in CR.prune.in > keep_sites.prune.in

# Prune from admixed individuals
plink2 --vcf $VCF \
--extract keep_sites.prune.in \
--keep $ALL_SAMPLES \
--recode vcf bgz
--allow-extra-chr \
--set-all-var-ids @:# \
--out ahmm_pruned_cbp

# convert VCF to Ancestry HMM input format
python3 vcf2ahmm.py -v ahmm_pruned_cbp.vcf.gz -s hmm_sample_mapping.txt