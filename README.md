# Internship_T2T

# Table of contents
1. [Introduction](#introduction)
2. [Usage](#Usage)
      - [Comparing SV indel balance](#SVcomparing)
      - [Comparing mismatch rate](#Mismatch)
      - [Comparing average coverage](#Coverage)


## Introduction <a name="introduction"></a>

### text prone to change
Over the years, several reference genomes for genomic analysis on humans have been developed with the most recent being the T2T-CHM13 (reference the paper) reference genome. Reference genomes are representative examples of the genes that a idealize organism of a certain species has. They are used to determine different kinds of variations based on the differences between the reference genome and the genomes of patients. Those deviations can then be used to diagnose genetic diseases or disorders, draw conclusions on the cause of those diseases and draw new genetic insights. As with any type of research the quality of the tools (including the reference genome) is important. Unfortunately, many difficult regions of the reference genome have been unresolved for decades now. This caused several genetic research studies to not include these regions, which might have influenced the results (reference the paper). However, do the differences between the new T2T-CHM13 (reference) and the currently used GRCh13 (reference) reference genomes actually warrant the T2T-CHM13 as a supplementary tool for the diagnostic analysis done by the Radboudumc?  


## Usage <a name="Usage"></a>

### Comparing SV indel balance <a name="SVcomparing"></a>

Follow these steps to make the lineplot showing the difference in SV indel balance between the GRCh38 and the T2T reference genomes:
1. Go to the server where the VCF files that you want to plot are located (one from the Hg38 reference genome and one from the T2T reference genome).
2. Run the [run_compare_indel_balance.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/run_compare_indel_balance.sh) bash script with first the input VCF file and then the output file. Do this for both the T2T and GRCh38 VCF files.
3. Run the [make_sv_lineplot.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/make_sv_lineplot.py) python script.
4. Give the path to the filtered T2T file.
5. Give the path to the filtered GRCh38 file.
6. The file SV_indel_comparison.jpg contains the line plot.

### Comparing the mismatch rate <a name="Mismatch"></a>

Follow these steps to make the scatterplot showing the difference in mismatch rate between the GRCh38 and the T2T reference genomes:
1. Get the Binary Alignment Map (BAM) files from the same sample where the reads are mapped against the GRCh38 or T2T reference genomes.
2. Run the [make_mismatch_plot.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/mismatch_plot/make_mismatch_plot.py) python script.
3. Give the path to the T2T BAM file.
4. Give the path to the GRCh38 BAM file.
5. The file mismatch_rate_comparison.jpg file contains the scatter plot.

### Comparing average coverage <a name="Coverage"></a>

Follow these steps to make the histogram showing the difference in average coverage calculated with regions consisting of 500 base pairs between the two reference genomes:
1. Get the same BAM files as used for the mismatch rate if using the same sample as used in that step is desired. Otherwise use BAM files from another sample that has BAM file from both reference genomes.
2. Run the [get_coverage_mean_regions.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/coverage_occurrence_histogram/get_coverage_mean_regions.sh) bash script with first the name of the prefix (name that will be given to sub output file), second the input BAM file and the name of the output file. Do this once with the T2T BAM file as input file and once with the GRCh38 as input file and change the prefix, so the same one won't be used twice.
3. Get the (name_prefix).regions.bed.gz file. This file contains the average coverage per region of 500 base pairs. One regions.bed.gz file for the T2T and GRCh38 genomes.
4. Run the [make_coverage_mean_histogram.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/coverage_occurrence_histogram/make_coverage_mean_histogram.py) python script.
5. Give the path to the T2T regions.bed.gz file.
6. Give the path to the GRCh38 regions.bed.gz file.
7. The file coverage_occurrences_histogram.jpg contains the histogram.
