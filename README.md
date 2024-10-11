# Internship_T2T

# Table of contents
1. [Introduction](#introduction)
2. [Usage](#Usage)
      - [Comparing SV indel balance](#SVcomparing)
      - [Comparing mismatch rate](#Mismatch)
      - [Comparing average coverage](#Coverage)
      - [Comparing low coverage regions](#Low_coverage)
      - [Comparing difficult regions with the low coverage regions](#)


## Introduction <a name="introduction"></a>

### text prone to change
Over the years, several reference genomes for genomic analysis on humans have been developed with the most recent being the T2T-CHM13 (reference the paper) reference genome. Reference genomes are representative examples of the genes that a idealize organism of a certain species has. They are used to determine different kinds of variations based on the differences between the reference genome and the genomes of patients. Those deviations can then be used to diagnose genetic diseases or disorders, draw conclusions on the cause of those diseases and draw new genetic insights. As with any type of research the quality of the tools (including the reference genome) is important. Unfortunately, many difficult regions of the reference genome have been unresolved for decades now. This caused several genetic research studies to not include these regions, which might have influenced the results (reference the paper). However, do the differences between the new T2T-CHM13 (reference) and the currently used GRCh13 (reference) reference genomes actually warrant the T2T-CHM13 as a supplementary tool for the diagnostic analysis done by the Radboudumc?  


## Usage <a name="Usage"></a>

### Comparing SV indel balance <a name="SVcomparing"></a>

Follow these steps to make the lineplot showing the difference in SV indel balance between the GRCh38 and the T2T reference genomes:
1. Go to the server where the VCF files that you want to plot are located (one from the Hg38 reference genome and one from the T2T reference genome).
2. Run the [get_amount_indels_from_file.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/sv_indel_balance_plot/get_amount_indels_from_file.sh) bash script with first the input VCF file and then the output file. Do this for both the T2T and GRCh38 VCF files. Following example is for the T2T VCF file:
```
bash get_amount_indels_from_file.sh T2T.vcf T2T_svs.txt
```
3. Run the [make_sv_lineplot.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/make_sv_lineplot.py) python script following the example:
```
python make_sv_lineplot.py t2t_sv.txt GRCh38_sv.txt
```
4. The file SV_indel_comparison.png contains the line plot.

### Comparing the mismatch rate <a name="Mismatch"></a>

Follow these steps to make the scatterplot showing the difference in mismatch rate between the GRCh38 and the T2T reference genomes:
1. Get the Binary Alignment Map (BAM) files from the same sample where the reads are mapped against the GRCh38 or T2T reference genomes.
2. Run the [make_mismatch_plot.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/mismatch_plot/make_mismatch_plot.py) python script following the example:
```
python make_coverage_plot.py T2T.bam GRCh38.bam
```
3. The file mismatch_rate_comparison.png file contains the scatter plot.

### Comparing average coverage <a name="Coverage"></a>

Follow these steps to make the histogram showing the difference in average coverage calculated with regions consisting of 500 base pairs between the two reference genomes:
1. Get the same BAM files as used for the mismatch rate if using the same sample as used in that step is desired. Otherwise use BAM files from another sample that has BAM file from both reference genomes.
2. Run the [get_coverage_mean_regions.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/coverage_occurrence_histogram/get_coverage_mean_regions.sh) bash script with first the name of the prefix (name that will be given to sub output file), second the input BAM file and the name of the output file. Do this once with the T2T BAM file as input file and once with the GRCh38 as input file and change the prefix, so the same one won't be used twice. Following example is for the T2T bam file:
```
python get_coverage_mean_regions.sh name_t2t_prefix T2T.bam
```
3. Get the (name_prefix).regions.bed.gz file. This file contains the average coverage per region of 500 base pairs. One regions.bed.gz file for the T2T and GRCh38 genomes.
4. Run the [make_coverage_mean_histogram.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/coverage_occurrence_histogram/make_coverage_mean_histogram.py) python script following the example:
```
python make_coverage_mean_histogram.py T2T.bed.gz GRCh38.bed.gz
```
5. The file coverage_occurrences_histogram.png contains the histogram.

### Comparing the regions with low coverage <a name="Low_coverage"></a>

Follow these steps to compare the regions with low coverage between the two reference genomes in the first 5 mega base pairs of chromosome 1:
1. Run the [get_low_coverage_regions.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/low_coverage_comparison/get_low_coverage_regions.sh) bash script with the two regions.bed.gz files obtained from the [get_coverage_mean_regions.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/coverage_occurrence_histogram/get_coverage_mean_regions.sh) with the following code:
```
bash get_low_coverage_regions.sh T2T_regions.bed.gz T2T_low_coverage_regions.bed GRCh38_regions.bed.gz GRCh38_low_coverage_regions.bed
```
2. Lift the coordinates of the low average regions in the just gotten output_T2T.bed over to the GRCh38 reference genome with your preferred liftover tool. This results in a bed file containing the coordinates for the low average coverage regions on the GRCh38 reference genome.
3. Run the [comparing_low_coverage_regions.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/low_coverage_comparison/comparing_low_coverage_regions.py) python script with the input being the T2T.bed, GRCh38.bed and T2T_to_GRCh38.bed.
```
python comparing_low_coverage_regions.py T2T_BED_file GRCh38_BED_file Liftover_BED_file
```
4. The file low_coverage_comparison.png contains the scatterplot.

### Comparing overlapping difficult regions with the low coverage regions <a name="category"></a>

Follow these steps to get the scatterplot in which the amount of times a low coverage region overlaps with a centromere, transcript, coding sequence (CDS) and/or structural duplication. These steps are written with the assumption that step 1 and 2 of comparing low coverage regions has already been done and the files have been made.
1. Obtain the hg38.ncbiRefSeq.gtf file [link](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/), a bed file with structural duplications and a gz file with the centromeres.
2. run the [make_low_coverage_categories_file.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/low_coverage_comparison/make_low_coverage_categories_file.sh) bash file like the following example:
```
bash make_low_coverage_categories_file.sh hg38.ncbiRefSeq.gtf hg38_segmental_duplications centromeres.txt.gz ouput_CDS_bed_file output_transcript_bed_file output_centromeres_bed_file T2T_low_coverage_regions.bed GRCh38_low_coverage regions.bed T2T_output.bed GRCh38_output.bed
```   
2. Open the [low_coverage regions.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/low_coverage_comparison/low_coverage_categories.py) python script, change the following line to the appropriate failed liftover value and save the change:
```
t2t_count_values_dict["Failed liftover"] = 207089
```
3. run the python like the following code example:
```
python script intersected_T2T_liftover.bed intersected_GRCh38.bed
```
bash /ifs/home/wout/PycharmProjects/Internship_T2T/low_coverage_comparison/make_low_coverage_categories_file.sh hg38.ncbiRefSeq.gtf hg38_segmental_duplications.bed centromeres.txt.gz test_bash_CDS test_bash_transcript test_bash_centromeres /ifs/data/research/projects/wout/projects/wp1_depth/data/BAMs_T2T_P3-D10/T2T_entire_genome_low_coverage.bed /ifs/data/research/projects/wout/projects/wp1_depth/data/BAMs_GRCh38_P3-D10/GRCh38_entire_genome_P3-D10_low_coverage.bed test_bash_T2T test_bash_GRCh38

