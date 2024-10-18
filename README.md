# Internship_T2T

# Table of contents
1. [Introduction](#introduction)
2. [Usage](#Usage)
      - [Comparing SV indel balance](#SVcomparing)
      - [Comparing mismatch rate](#Mismatch)
      - [Comparing average coverage](#Coverage)
      - [Comparing low coverage regions](#Low_coverage)
      - [Comparing difficult regions with the low coverage regions](#category)


## Introduction <a name="introduction"></a>

### text prone to change
Over the years, several reference genomes for genomic analysis on humans have been developed with the most recent being the T2T-CHM13 (reference the paper) reference genome. Reference genomes are representative examples of the genes that a idealize organism of a certain species has. They are used to determine different kinds of variations based on the differences between the reference genome and the genomes of patients. Those deviations can then be used to diagnose genetic diseases or disorders, draw conclusions on the cause of those diseases and draw new genetic insights. As with any type of research the quality of the tools (including the reference genome) is important. Unfortunately, many difficult regions of the reference genome have been unresolved for decades now. This caused several genetic research studies to not include these regions, which might have influenced the results (reference the paper). However, do the differences between the new T2T-CHM13 (reference) and the currently used GRCh13 (reference) reference genomes actually warrant the T2T-CHM13 as a supplementary tool for the diagnostic analysis done by the Radboudumc?  


## Usage <a name="Usage"></a>

### Comparing SV indel balance <a name="SVcomparing"></a>

Follow these steps to create a line plot comparing the indel balance between the GRCh38 and T2T reference genomes, and a bar plot showing the difference in the total number of insertions and deletions between the two genomes:
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
5. The file Total_indel_comparison.png contains the barplot.

### Comparing the mismatch rate <a name="Mismatch"></a>

Follow these steps to make the scatterplot showing the difference in mismatch rate between the GRCh38 and the T2T reference genomes:
1. Get the Binary Alignment Map (BAM) files from the same sample where the reads are mapped against the GRCh38 or T2T reference genomes.
2. Run the [make_mismatch_plot.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/mismatch_plot/make_mismatch_plot.py) python script following the example:
```
python make_coverage_plot.py T2T.bam GRCh38.bam
```
3. The file mismatch_rate_comparison.png file contains the scatter plot.

### Comparing average coverage <a name="Coverage"></a>

Follow these steps to create a histogram comparing the difference in average coverage for 500 base-pair regions between two reference genomes, and a scatterplot showing the difference in standard deviation of average coverage between the same regions:
1. Get the same BAM files as used for the mismatch rate if using the same sample as used in that step is desired. Otherwise use BAM files from another sample that has BAM file from both reference genomes.
2. Run the [get_coverage_mean_regions.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/coverage_occurrence_histogram/get_coverage_mean_regions.sh) bash script with first the name of the prefix (name that will be given to sub output file), second the input BAM file and the name of the output file. Do this once with the T2T BAM file as input file and once with the GRCh38 as input file and change the prefix, so the same one won't be used twice. Following example is for the T2T bam file:
```
python get_coverage_mean_regions.sh T2T.bam name_prefix
```
3. Get the (name_prefix).regions.bed.gz file. This file contains the average coverage per region of 500 base pairs. One regions.bed.gz file for the T2T and GRCh38 genomes.
4. Run the [make_coverage_mean_histogram.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/coverage_occurrence_histogram/make_coverage_mean_histogram.py) python script following the example:
```
python make_coverage_mean_histogram.py T2T.bed.gz GRCh38.bed.gz
```
5. The file coverage_occurrences_histogram.png contains the histogram.
6. The file standard_deviation_coverage.png contains the scatterplot

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

Follow these steps to create a bar plot showing the frequency of low coverage regions overlapping with centromeres, transcripts, coding sequences (CDS), and/or structural duplications. The steps for T2T and GRCh38 are different because they require different input files and are thus written separately. These steps are written with the assumption that step 1 of comparing low coverage regions has already been done and the low coverage files have been made.

First the GRCh38:
1. Obtain the hg38.ncbiRefSeq.gtf file ([link](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/)), a bed file with structural duplications (file was given directly, so no example is available) and a gz file with the centromeres ([link](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/)).
2. Run the [make_low_coverage_categories_hg38.sh](htts://github.com/WoutPoelen/Internship_T2T/blob/main/low_coverage_comparison/make_low_coverage_categories_hg38.sh) bash file like the following example:
```
bash make_low_coverage_categories_hg38.sh hg38.ncbiRefSeq.gtf hg38_segmental_duplications.bed centromeres.txt.gz ouput_CDS.bed output_transcript.bed output_centromeres.bed GRCh38_low_coverage.regions.bed GRCh38_output.bed
```

Then the T2T-CHM13:
1. Obtain the refseq t2t genes gtf file ([link](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/)).
2. run [convert_T2T_genes_file.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/low_coverage_comparison/convert_T2T_genes_file.sh) to get the normal chromosome names.
3. Obtain the bigbed (bb) file containing pericentromeric and centromeric satellites ([link](https://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=345820279_xEDUaM4aXhxuQpQp1EiinRxuQAFH&db=hub_567047_hs1&c=chr9&g=hub_567047_censat)). Warning!!! This file doesn't contain every centromere perfectely. Some chromosomes have regions in their centromeres that will not be compared. This problem will be resolved once there is a file which contains all of the regions.
4. Change the bb file into a bed file with the following code:
```bigBedToBed censat.bb name_file.bed```
5. Obtain a bed file with the structural duplications (file was given directly, so no example is available).
6. Run the [make_low_coverage_categories_t2t.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/low_coverage_comparison/make_low_coverage_categories_t2t.sh) bash file like the following example:
```
bash make_low_coverage_categories_t2t.sh converted_refseq_t2t.gtf Structural_Duplications.bed centromeres_CenSat_Annotation.bed CDS_output.bed Transcript_output.bed Centromeres_output.bed T2T_low_coverage.bed T2T_output_categorical.bed
```

To make the plot:
1. Run the [low_coverage regions.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/low_coverage_comparison/low_coverage_categories.py) python script like the following code example:
```
python script T2T_output_categorical.bed GRCh38_output_categorical.bed
```
2. The barplot is saves as low_coverage_categories_barplot.png

