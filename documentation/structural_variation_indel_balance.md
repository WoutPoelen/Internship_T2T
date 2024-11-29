### Calculating and comparing the SV indel balance

#### Introduction
This step is performed to see if the indel balance is better in the T2T-CHM13 reference genome. 
This shows the reliability of the new reference genome when used for structural variant calling compared
to the GRCh38 reference genome. Structural variation in the genome can distort the read mapping in and around
a gene causing genetic differences between the sample and the reference genome in the region to become untrustworthy. 
This can have a great impact on further steps of genetic analysis and diagnosis.

#### Steps to perform

Follow these steps to create a line plot comparing the indel balance between the GRCh38 and T2T reference genomes, and a bar plot showing the difference in the total number of insertions and deletions between the two genomes:
1. Go to the server where the VCF files that you want to plot are located (one from the Hg38 reference genome and one from the T2T reference genome).
2. Run the [get_amount_indels_from_file.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/bash/get_amount_indels_from_file.sh) bash script with first the input VCF file and then the output file. Do this for both the T2T and GRCh38 VCF files. Following example is for the T2T VCF file:
```
bash get_amount_indels_from_file.sh T2T.vcf T2T_svs.txt
```
3. Run the [structural_variation_lineplot.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/python/structural_variation_lineplot.py) python script following the example:
```
python make_sv_lineplot.py t2t_sv.txt GRCh38_sv.txt
```
4. The file SV_indel_comparison.png contains the line plot.
5. The file Total_indel_comparison.png contains the barplot.

#### Advised next step
Getting the mean mapping coverage and plotting it into a histogram 
[mapping_coverage_mean_histogram.md](https://github.com/WoutPoelen/Internship_T2T/tree/main/documentation/mapping_coverage_mean_histogram.md)