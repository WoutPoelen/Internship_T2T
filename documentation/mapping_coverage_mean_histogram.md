### Comparing the average coverage of 500bp regions 

---
#### Introduction
This step is performed to create a histogram comparing the difference in average coverage for 500 base-pair regions between 
two reference genomes, and a scatterplot showing the difference in standard deviation of average coverage between the same regions.

---
#### Steps to perform
1. Get the same BAM files as used for the mismatch rate if using the same sample as used in that step is desired. Otherwise, use BAM files from another sample that has BAM file from both reference genomes.
2. Run the [get_coverage_mean_regions.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/bash/get_coverage_mean_regions.sh) 
bash script with first the name of the prefix (name that will be given to sub output file), second the input BAM file 
and the name of the output file. Do this once with the T2T BAM file as input file and once with the GRCh38 as input file 
and change the prefix, so the same one won't be used twice. Following example is for the T2T bam file:
```
bash get_coverage_mean_regions.sh T2T.bam name_prefix
```
3. Get the (name_prefix).regions.bed.gz file. This file contains the average coverage per region of 500 base pairs. 
One regions.bed.gz file for the T2T and GRCh38 genomes.
4. Run the [mapping_coverage_mean_histogram.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/python/mapping_coverage_mean_histogram.py) python script following the example:
```
python make_coverage_mean_histogram.py T2T.bed.gz GRCh38.bed.gz
```
5. The file coverage_occurrences_histogram.png contains the histogram.
This histogram also contains the mean and standard deviation of the coverage. The occurrences of regions with coverage 
higher than 80 are combined and plotted at the 80+ label.
6. The file standard_deviation_coverage.png contains the scatterplot.

---
#### Advised next step
Plotting the 500bp regions whose read mapping coverage is below a to be calculated threshold in karyoplots
[low_coverage_karyoplot.md](https://github.com/WoutPoelen/Internship_T2T/tree/main/documentation/low_coverage_karyoplot.md)

or

Separating the low coverage regions into four categories depending on their location (Centromeres, coding sequences, 
transcripts and segmental duplications) [separating-low_coverage_into_categories.md](https://github.com/WoutPoelen/Internship_T2T/tree/main/documentation/separating_low_coverage_into_categories.md).