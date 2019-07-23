#Co-Cleaning Command Line Parameters - Uses tumor and normal bam files
#Step one - RealignTargetCreator
java -jar "/home/users/industry/chicago-university/meghan/miniconda3/bin/GenomeAnalysisTK.jar" -T RealignerTargetCreator -R  "/home/users/industry/chicago-university/meghan/GDC-Pipeline/GRCh38.d1.vd1.fa" -known "/home/users/industry/chicago-university/meghan/GDC-Pipeline/common_all_20180418.vcf" -I "/home/users/industry/chicago-university/meghan/GDC-Pipeline/S5-OUTPUT/Sample1/Tumor/S5S1Tumor.bam" -I "/home/users/industry/chicago-university/meghan/GDC-Pipeline/S5-OUTPUT/Sample1/Normal/S5S1Norm.bam" -o "/home/users/industry/chicago-university/meghan/GDC-Pipeline/STCC-OUTPUT/S1Output.intervals"

#Step two - IndelRealigner
#java -jar "/home/users/industry/chicago-university/meghan/miniconda3/bin/GenomeAnalysisTK.jar" -T IndelRealigner -R "/home/users/industry/chicago-university/meghan/GDC-Pipeline/GRCh38.d1.vd1.fa" -known "/home/users/industry/chicago-university/meghan/GDC-Pipeline/common_all_20180418.vcf" -targetIntervals "/home/users/industry/chicago-university/meghan/GDC-Pipeline/STCC-OUTPUT/S1Output.intervals" --noOriginalAlignmentTags -I "/home/users/industry/chicago-university/meghan/GDC-Pipeline/S5-OUTPUT/Sample1/Tumor/S5S1Tumor.bam" -I "/home/users/industry/chicago-university/meghan/GDC-Pipeline/S5-OUTPUT/Sample1/Normal/S5S1Norm.bam" -o "/home/users/industry/chicago-university/meghan/GDC-Pipeline/STCC-OUTPUT/output.map"

#Step three - BaseRecalibrator; dbSNP v.144
#java -jar "/home/users/industry/chicago-university/meghan/miniconda3/bin/GenomeAnalysisTK.jar" -T BaseRecalibrator -R "/home/users/industry/chicago-university/meghan/GDC-Pipeline/GRCh38.d1.vd1.fa" -I "/home/users/industry/chicago-university/meghan/GDC-Pipeline/S5-OUTPUT/Sample1/Tumor/S5S1Tumor.bam" -I "/home/users/industry/chicago-university/meghan/GDC-Pipeline/S5-OUTPUT/Sample1/Normal/S5S1Norm.bam" -knownSites "/home/users/industry/chicago-university/meghan/GDC-Pipeline/common_all_20180418.vcf" -o "/home/users/industry/chicago-university/meghan/GDC-Pipeline/STCC-OUTPUT/S1bqsr.grp"

#Step four - PrintReads
#java -jar "/home/users/industry/chicago-university/meghan/miniconda3/bin/GenomeAnalysisTK.jar" -T PrintReads -R "/home/users/industry/chicago-university/meghan/GDC-Pipeline/GRCh38.d1.vd1.fa" -I "/home/users/industry/chicago-university/meghan/GDC-Pipeline/S5-OUTPUT/Sample1/Tumor/S5S1Tumor.bam" -I "/home/users/industry/chicago-university/meghan/GDC-Pipeline/S5-OUTPUT/Sample1/Normal/S5S1Norm.bam" --BQSR "/home/users/industry/chicago-university/meghan/GDC-Pipeline/STCC-OUTPUT/S1bqsr.grp" -o "/home/users/industry/chicago-university/meghan/GDC-Pipeline/STCC-OUTPUT/S1Output.bam"

