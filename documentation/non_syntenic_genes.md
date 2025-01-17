### Comparing genes with non-syntenic regions

---

#### Introduction
This step is performed to get the genes which are located in non-syntenic regions in T2T
compared to GRCh38. These genes contain sequences that aren't known in GRCh38 and could thus be interesting
if medically relevant because they can then be used for further research and genetic diagnostics. 

---

#### Steps to perform
**269 genes**
1. Retrieve the 269 medically relevant genes from the [T2T article](https://www.science.org/doi/10.1126/science.abl3533#supplementary-materials) 
by downloading the supplementary tables and selecting only the genes in Table S13 that have an SNV benchmark in T2T
2. Remove the genes that don't have T2T coordinates by running the following command:
```
grep -v '^\s' medically_relevant_genes.bed | awk '{OFS]"\t";$1=$1}' > removed_empty_coordinates_genes.bed
```
3. Overlap the remaining medically relevant genes with the non-syntenic regions by running the following command:
```
bedtools intersect -a removed_empty_coordinates_genes.bed -b non-syntenic_regions.bed > output_file
```

**All genes**
1. Overlap the coding sequences of T2T (obtained in [make_low_coverage_categories_t2t.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/bash/make_low_coverage_categories_t2t.sh))
by running the following command:
```
bedtools intersect -a T2T_CDS_regions -b non_syntenic_regions.bed | sort | uniq > output_file
```
2. Overlap the output file with the segmental duplications to get the genes whose results are less trustworthy 
by running the following command:
```
bedtools intersect -a all_genes_intersected.bed -b segmental_duplications.bed | cut -f 4 | uniq > all_CDS_overlapping_non_syntenic_SD.bed
```
This results in a list with genes that overlap with a non-syntenic and segmental duplicated regions. These genes are not
really trustworthy enough to draw conclusions from.

3. Get the genes that don't overlap with a SD by overlapping the step 1 output file with the segmental like the 
following command:
```
bedtools intersect -v -a all_genes_intersected.bed -b segmental_duplications.bed | cut -f 4 | uniq > all_CDS_overlapping_non_syntenic_SD.bed
```
This results in a list with genes that overlap with a non-syntenic but not with a segmental duplicated regions. 
These genes are more trustworthy to draw conclusions from and could potentially be medically relevant.

---

#### Possible next steps
Looking at the functions of the genes that are not overlapping with a segmental duplications, but are in a non-syntenic 
regions to see if they are medically relevant genes.