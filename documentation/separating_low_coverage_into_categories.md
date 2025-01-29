### Overlapping regions that are difficult to map and regions with low coverage

---
#### Introduction
This step is performed to create a bar plot showing difficult to sequence regions overlap with the low coverage 
regions. Specifically looking at centromeres, transcripts, coding sequences (CDS), and/or structural duplications. 
This is done to see the read mapping coverage improved in the coding sequences. Since these regions are
important to see if there are genes which improved in the new reference genome. The other three categories are used
because the also have an impact on the genes which reside in that region.

---
#### Steps to perform

**Obtaining the low coverage regions:**

Run the [obtaining_low_coverage_regions.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/python/obtaining_low_coverage_regions.py) 
python script with the two regions.bed.gz files obtained from the
[get_coverage_mean_regions.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/bash/get_coverage_mean_regions.sh) with the following code:

```
python obtaining_low_coverage_regions.py T2T_regions.bed.gz GRCh38_regions.bed.gz T2T_low_coverage_regions.bed GRCh38_low_coverage_regions.bed
```

This script calculates the median for the autosomal chromosomes and the sex chromosomes separatly (alternative sequences were not used). 
Then calculates the low coverage threshold by dividing the median for the autosomal chromosomes by three and the median for the sex chromosomes by six.


**Identifying low coverage regions in GRCh38 that overlap with difficult-to-sequence areas:**
1. Obtain the hg38.ncbiRefSeq.gtf file ([link](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/)), a bed 
file with structural duplications (file was given directly, so no example is available) and a gz file with the centromeres ([link](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/)).
2. Run the [make_low_coverage_categories_hg38.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/bash/make_low_coverage_categories_hg38.sh) bash file like the following example:
```
bash make_low_coverage_categories_hg38.sh hg38.ncbiRefSeq.gtf hg38_segmental_duplications.bed centromeres.txt.gz ouput_CDS.bed output_transcript.bed output_centromeres.bed GRCh38_low_coverage.regions.bed GRCh38_output.bed
```

**Identifying low coverage regions in T2T that overlap with difficult-to-sequence areas:**

1. Obtain the refseq t2t genes gtf file [link](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/).
2. run [convert_T2T_genes_RefSeq.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/bash/convert_T2T_genes_RefSeq.sh) to get the normal chromosome names.
3. Obtain the bigbed (bb) file containing peri-centromere and centromere satellites ([link](https://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=345820279_xEDUaM4aXhxuQpQp1EiinRxuQAFH&db=hub_567047_hs1&c=chr9&g=hub_567047_censat)). 
Warning!!! This file doesn't contain every centromere perfectly. Some chromosomes have regions in their centromeres that will not be compared. This problem will be resolved once there is a file which contains all the regions.
4. Change the bb file into a bed file with the following code:
```bigBedToBed censat.bb name_file.bed```
5. Obtain a bed file with the structural duplications (file was given directly, so no example is available).
6. Run the [make_low_coverage_categories_t2t.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/bash/make_low_coverage_categories_t2t.sh) bash file like the following example:
```
bash make_low_coverage_categories_t2t.sh converted_refseq_t2t.gtf Structural_Duplications.bed centromeres_CenSat_Annotation.bed CDS_output.bed Transcript_output.bed Centromeres_output.bed T2T_low_coverage.bed T2T_output_categorical.bed
```

**To make the plot:**
1. Run the [separating_low_coverage_into_categories.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/python/separating_low_coverage_into_categories.py) python script like the following code example:
```
python script T2T_output_categorical.bed GRCh38_output_categorical.bed
```
2. The barplot is saves as low_coverage_categories_barplot.png.

---
#### Advised next step
Getting the genes who has a coding sequence overlap with a low coverage region [comparing_low_coverage_coding_sequence.md](https://github.com/WoutPoelen/Internship_T2T/tree/main/documentation/comparing_low_coverage_coding_sequence.md)
