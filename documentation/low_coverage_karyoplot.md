### Low coverage comparison with a karyoplot 

#### Introduction
This step is performed to make karyoplots of the low coverage regions in the GRCh38 and T2T reference genome.
This is done to get a clearer view of the locations and the amount of low coverage regions in the reference genomes. 

#### Steps to perform

**Obtaining the low coverage regions:**

Run the [obtaining_low_coverage_regions.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/python/obtaining_low_coverage_regions.py) 
python script with the two regions.bed.gz files obtained from the
[get_coverage_mean_regions.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/bash/get_coverage_mean_regions.sh) with the following code:

The Rscripts cannot be called on the command line.

**Installing necessary packages:**
1. Run the following code to install BiocManager
```
Install.packages("BiocManager")
```
2. Run the following code to install KaryoploteR
```
BiocManager::install("KaryoploteR")
```

**GRCh38 karyoplot:**
1. Copy the [karyoplot_GRCh38.R](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/r/karyoplot_GRCh38.R) R script and paste it in Rstudio.
2. Import the low coverage GRCh38 bed file (made in [obtain_low_coverage_regions.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/python/obtaining_low_coverage_regions.py))  as a dataset. This will automatically make it a table which is necessary for the script to work.
3. Change the following lines to the name of the GRCh38 low coverage file:
```
(chr=GRCh38_low_coverage_HG002["V1"],
start=GRCh38_low_coverage_HG002["V2"], 
end=GRCh38_low_coverage_HG002["V3"],
y=GRCh38_low_coverage_HG002["V4"]))
```
4. Run the script.
5. Export the karyoplot to a specific location.


**T2T karyoplot:**

The T2T reference genome isn't available in KaryoplotR (version 1.30.0). So a custom reference genome needs to be created.
1. Copy the [karyoplot_T2T.R](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/r/karyoplot_T2T.R) R script and paste it in Rstudio
2. Import the low coverage T2T bed file (made in [obtain_low_coverage_regions.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/python/obtaining_low_coverage_regions.py)) as a dataset. This will automatically make it a table which is necessary for the script to work.
3. Change the following lines to the name of the T2T low coverage file:
```
chr=T2T_low_coverage_HG002["V1"], 
start=T2T_low_coverage_HG002["V2"], 
end=T2T_low_coverage_HG002["V3"],
y=T2T_low_coverage_HG002["V4"])
```
4. Get the chromosomes length by downloading the [NIH T2T chromosomes dataset](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/) at the bottom of the website and importing the file as a dataset.
5. Change the following lines to the name of the newly imported chromosome file (start doesn't have to change):
```
chr=chromosome_report_t2t$UCSC.style.name,
start=0, 
end=chromosome_report_t2t$Seq.length
```
6. Get the centromere locations by importing the bed file made in step 4 of the T2T part of [Comparing difficult regions with the low coverage regions](#category)
7. Change the following lines to the name of the Censat centromere file:
```
chr=filtered_centromeres_CenSat$V1,
start=filtered_centromeres_CenSat$V2,
end=filtered_centromeres_CenSat$V3)
```
8. Run the script.
9. Export the karyoplot to a specific location.

#### Advised next step
Separating the low coverage regions into four categories depending on their location
(Centromeres, coding sequences, transcripts and segmental duplications):
[separating-low_coverage_into_categories.md](https://github.com/WoutPoelen/Internship_T2T/tree/main/documentation/separating_low_coverage_into_categories.md)
