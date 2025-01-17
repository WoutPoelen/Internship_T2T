### Looking at the non-syntenic structural variants from multiple samples

---

#### Introduction
This step is performed to identify the structural variants present in the T2T reference genome that are absent in 
GRCh38 and are overlapping with a coding sequence. These structural variants can be used as an argument towards using T2T instead of GRCh38 for genomic diagnosis,
because they could provide an explanation of a previously unknown reason for a genes loss of function.

---

#### Steps to perform
1. Run the [T2T_multiple_sample_vcf_files_processing.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/bash/T2T_multiple_sample_vcf_files_processing.sh)
bash script like the following example:
```
bash T2T_multiple_sample_vcf_files_processing.sh non_syntenic_regions_file coding_sequences_T2T_file path_directory_samples path_output_directory_bcftools path_output_gene_intersection path_output_SV_intersection
```
This results in three directories being filled:

- path_output_directory_bcftools contains a separate file for each sample containing the insertions and deletions from that sample.
- path_output_gene_intersection contains a separate file for each sample containing the coding sequences (and gene) that overlap with an indel and a non-syntenic region. 
- path_output_SV_intersection contains a separate file for each sample containing the indels in that sample that overlap with a non-syntenic region and a coding sequence. 
Each indel gets written once even if it overlaps with multiple CDS.

2. Run the following command in the output_SV_intersection directory to get a list of the amount of SVs in non-syntenic and coding sequences:
```
wc -l * | sort -nr | tr -s " " > SVs_samples_occurences.txt
```

3. Go through the intersected SV files created in step 1 and make separate files (one for deletions and one for insertions) containing the amount of insertions or deletions each sample has to get a better view over a possible bias.

Deletions
```
for file in *.bed; do 
    filename=$(basename "$file" .bed)
    echo -n "$filename"
    grep -o -i 'DEL' "$file" | wc -l
done > sort -nr > DEL_output_file_name.txt
```
Insertions
```
for file in *.bed; do 
    filename=$(basename "$file" .bed)
    echo -n "$filename"
    grep -o -i 'DEL' "$file" | wc -l
done > sort -nr > INS_output_file_name.txt
```

4. Run the [multiple_sample_non_syntenic_SV.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/python/multiple_sample_non_syntenic_SV.py)
python script with the following code:
```
python multiple_sample_non_syntnenic_SV.py DEL_output_file_name.txt INS_output_file_name.txt
```
The barplot is saved as: Multiple_sample_non_syntenic_SV_occurences.png


---

#### Possible next steps
1. Count in how many files a gene overlaps with a SV and non-syntenic regions in T2T using the following code example:
```
for file in path_output_gene_intersection/*; do awk '{print $4}' "$file" | sort -u; done | sort | uniq -c | sort -nr > outputfile
```
2. Run the [GRCh38_multiple_sample_vcf_files_processing.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/bash/GRCh38_multiple_sample_vcf_files_processing.sh)
bash script like the following example to get the SVs in GRCh38 that overlap with a coding sequence.
```
bash GRCh38_multiple_sample_vcf_files_processing.sh non_syntenic_regions_file coding_sequences_T2T_file path_directory_samples path_output_directory_bcftools path_output_bedtools
```
This results in two directories being filled:

- path_output_directory_bcftools contains a file for each sample containing the insertions and deletions from that sample.
- path_output_directory_bedtools contains a file for each sample containing the indels that overlap with a coding sequence.
