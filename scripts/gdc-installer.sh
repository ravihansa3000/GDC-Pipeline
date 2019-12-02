#!/bin/bash
set -e
############################################
# Install script for GDC Pipeline Workflow #
############################################

CONDA_ENV_GDC='gdc'
OUTPUT_DIR=${1:-~/Downloads}
mkdir -p $OUTPUT_DIR

conda update -y -n root conda
ENVS=$(conda env list | awk '{print $1}')
source ~/anaconda3/etc/profile.d/conda.sh
if [[ $ENVS = *"$CONDA_ENV_GDC"* ]]; then
    conda activate $CONDA_ENV_GDC
else 
    conda create -y --name $CONDA_ENV_GDC python=3.7
    conda activate $CONDA_ENV_GDC
fi;

# Install dependencies
conda install pip
pip install setuptools

# Parsl framework
# pip install parsl==0.9.0

# Install latest Parsl - required for Python 3.7+
[ -d "$OUTPUT_DIR/parsl" ] && rm -rf "$OUTPUT_DIR/parsl"
git clone https://github.com/Parsl/parsl.git "$OUTPUT_DIR/parsl"
cd $OUTPUT_DIR/parsl
python setup.py install
cd -

# SWAG framework
[ -d "$OUTPUT_DIR/Swag" ] && rm -rf "$OUTPUT_DIR/Swag"
git clone https://github.com/PittGenomics/Swag.git "$OUTPUT_DIR/Swag"
cd $OUTPUT_DIR/Swag
python setup.py install
cd -

# Parsl dependencies
pip install --no-binary pyzmq pyzmq
pip install pysqlite3
pip install sqlalchemy
pip install sqlalchemy_utils

### biobambam ver. 2.0.87 -updated to most recent version
conda install -y -c bioconda biobambam=2.0.87

### bwa ver. 0.7.15
conda install -y -c bioconda  bwa=0.7.15

### samtools ver. 1.1
# GDC Pipeline VARSCAN: Step 1: Mpileup; Samtools 1.1
conda install -y -c bioconda samtools=1.1
cp $CONDA_PREFIX/bin/samtools $CONDA_PREFIX/bin/samtools-1.1


### samtools ver. 1.3.1
# GDC Pipeline DNA Seq alignment: STEP 2: BWA ALIGNMENT - BWA 0.7.15 - SAMTOOLS 1.3.1
conda install -y -c bioconda samtools=1.3.1
cp $CONDA_PREFIX/bin/samtools $CONDA_PREFIX/bin/samtools-1.3.1

### picard ver. 2.5.0
conda install -y -c bioconda picard=2.5.0

### somatic-sniper ver. 1.0.5.0
conda install -y -c bioconda somatic-sniper=1.0.5.0

### varscan ver. 2.4.0
conda install -y -c bioconda varscan=2.4.0

### GATK ver. 4.0.4.0
# only used in tumor-only workflow
conda install -y -c bioconda gatk4=4.0.4.0

### MuSE ver. 1.0
conda install -y -c bioconda muse=1.0.rc

### GATK ver. 3.5-5 GDC ver. nightly-2016-02-25gf39d340
#FIXME: does not give GenomeAnalysis.tk file due to licensing
conda install -y -c bioconda gatk=3.5

wget -O $OUTPUT_DIR/GATK3-3.5_GenomeAnalysisTK.tar.bz2 "https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.5-0-g36282e4"
$CONDA_PREFIX/bin/../opt/gatk-3.5/gatk-register.sh $OUTPUT_DIR/GATK3-3.5_GenomeAnalysisTK.tar.bz2

### Variant Annotation Workflow (VEP)
conda install -y -c bioconda ensembl-vep=88

printf "=== success ==="
