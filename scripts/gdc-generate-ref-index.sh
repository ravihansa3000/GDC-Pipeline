#!/bin/bash
set -e
###############################################################
# Generate genome reference indexes for GDC Pipeline Workflow #
###############################################################

CONDA_ENV_GDC='gdc'
GENOME_REF_DIR=${1-'/media/akila/DATA/GDC/GDC-Pipeline-Data/reference'}
GENOME_REF="$GENOME_REF_DIR/GRCh38.d1.vd1.fa"

ENVS=$(conda env list | awk '{print $1}')
source ~/anaconda3/etc/profile.d/conda.sh
if [[ $ENVS = *"$CONDA_ENV_GDC"* ]]; then
    conda activate $CONDA_ENV_GDC
else 
    echo "creating new conda envrionment: $CONDA_ENV_GDC"
    conda create --name $CONDA_ENV_GDC
    conda activate $CONDA_ENV_GDC
fi;

# Install dependencies
conda install -y -c bioconda bcftools tabix vcftools
conda install -y -c bioconda samtools=1.3.1 picard=2.5.0

# creating index (fai)
printf "creating reference index for: $GENOME_REF \n"
samtools faidx $GENOME_REF

# creating dictionary (dict)
printf "creating dictionary for: $GENOME_REF \n"
PICARD_JAR="$CONDA_PREFIX/share/picard-2.5.0-2/picard.jar"
java -jar $PICARD_JAR CreateSequenceDictionary R=$GENOME_REF O="$GENOME_REF_DIR/GRCh38.d1.vd1.fa.dict"

# Combining Cosmic vcf files
# vcf-concat CosmicCodingMuts.vcf.gz CosmicNonCodingVariants.vcf.gz | gzip -c  > CosmicCoding.vcf.gz
# vcf-sort > CosmicCoding.vcf.gz

printf "=== success ==="