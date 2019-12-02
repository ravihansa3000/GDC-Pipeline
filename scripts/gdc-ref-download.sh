#!/bin/bash
set -e
#############################################
# Download script for GDC Pipeline Workflow #
#############################################

OUTPUT_DIR=${1:-'~/Downloads'}

### dbSNP ver. 144
wget -O "${OUTPUT_DIR}/dbsnp_144.hg38.vcf.gz" "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_144.hg38.vcf.gz"
wget -O "${OUTPUT_DIR}/dbsnp_144.hg38.vcf.gz.tbi" "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_144.hg38.vcf.gz.tbi"

### GDC co-cleaning workflow - known indels
wget -O "${OUTPUT_DIR}/Homo_sapiens_assembly38.known_indels.vcf.gz" "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz"
wget -O "${OUTPUT_DIR}/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi" "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"

### Cosmic.vcf - Catalogue Of Somatic Mutations In Cancer
### Make account here: https://cancer.sanger.ac.uk/cosmic
### Follow scripted directions to download vcf files for coding and non-coding mutations
### https://cancer.sanger.ac.uk/cosmic/file_download_info?data=GRCh38%2Fcosmic%2Fv89%2FVCF%2FCosmicCodingMuts.vcf.gz
### https://cancer.sanger.ac.uk/cosmic/file_download_info?data=GRCh38%2Fcosmic%2Fv89%2FVCF%2FCosmicNonCodingVariants.vcf.gz

### reference genome - GRCh38.d1.vd1
wget -O "${OUTPUT_DIR}/GRCh38.d1.vd1.fa.tar.gz" "https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834"
tar -zxf "${OUTPUT_DIR}/GRCh38.d1.vd1.fa.tar.gz" -C $OUTPUT_DIR

### indexed reference genome - GRCh38.d1.vd1_BWA
wget -O "${OUTPUT_DIR}/GRCh38.d1.vd1_BWA.tar.gz" "https://api.gdc.cancer.gov/data/25217ec9-af07-4a17-8db9-101271ee7225"
tar -zxf GRCh38.d1.vd1_BWA.tar.gz -C $OUTPUT_DIR

wget -O "${OUTPUT_DIR}/strelka-2.9.10.tar.gz" "https://github.com/Illumina/strelka/archive/v2.9.10.tar.gz"
tar -zxf "${OUTPUT_DIR}/strelka-2.9.10.tar.gz"

printf "=== success ==="