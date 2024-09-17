# Internship_T2T

# Table of contents
1. [Introduction](#introduction)
2. [Usage](#Usage)
      - [Comparing SV indel balance](#SVcomparing)


## Introduction <a name="introduction"></a>

### text prone to change
Over the years several reference genomes for humans have been developed with the most recent being the T2T-CHM13 reference genome. Reference genomes are representative examples of the genes that a idealize organism of a certain species has. They are used to determine different kinds of variations based on the differences between the reference genome and the genomes of patients. Those deviations can then be used to diagnose genetic diseases or disorders, draw conclusions on the cause of those diseases and draw new genetic insights. As with any type of research the quality of the tools (including the reference genome) is important. Unfortunately, many difficult regions of the reference genome have been unresolved for decades now. This caused several genetic research studies to not include these regions, which might have influenced the results. However, are the differences between the new T2T-CHM13 and the currently used GRCh13 reference genomes actually big enough to warrant a change in reference genome for the diagnostic analysis done by the Radboudumc?



## Usage <a name="Usage"></a>

### Comparing SV indel balance <a name="SVcomparing"></a>

Follow these steps to make the lineplot showing the difference in SV indel balance between the Hg38 and the T2T reference genomes:
1. Go to the server where the VCF files that you want to plot are located (one from the Hg38 reference genome and one from the T2T reference genome).
2. Perform the following bcftools command (assuming that bcftools is installed):
   
- bcftools query -f "%CHROM %POS %FILTER [%GT] %INFO/SVLEN\n" -i 'INFO/SVTYPE="DEL" || INFO/SVTYPE="INS"' Input_file.vcf > ~/output_file.txt

5. run the [make_sv_lineplot.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/make_sv_lineplot.py) script.
6. Give the path to the filtered T2T file
7. Give the path to the filtered Hg38 file

