### Genes of which the coding sequences are in low coverage regions 

---
#### Introduction
This step is performed to get the amount of genes whose coding sequence overlap with a low coverage region.
This is done to see if there are genes whose read mapping improved in these sequence. Since they are the part of the 
gene which actually codes for the protein and where a variant can thus have a significant impact on the function of the 
gene.  

---
#### Steps to perform
1. Run the [low_coverage_coding_sequence_regions.sh](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/bash/low_coverage_coding_sequence_regions.sh) bash script with these files in the same order(some files get written to during the execution of the script and are given as empty files):


| **File Name**                               | **Description**                                                                                                                                                                                                                                                                                               |
|---------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| T2T categorical file                        | BED file containing the number of difficult-to-map regions that overlap with certain low coverage regions in the T2T reference genome.                                                                                                                              |
| T2T categorical CDS only file               | Initially empty BED file, populated with low coverage regions overlapping one or more CDS regions in the GRCh38 reference genome during script execution.                                                                                                            |
| GRCh38 categorical file                     | BED file containing the number of difficult-to-map regions that overlap with certain low coverage regions in the GRCh38 reference genome.                                                                                                                           |
| GRCh38 categorical CDS only file            | Initially empty BED file, populated with low coverage regions overlapping one or more CDS regions in the T2T reference genome during script execution.                                                                                                               |
| T2T RefSeq file                             | GTF file containing the genes in the T2T reference genome, obtained from step 2 of the T2T intersection steps in the linked section.                                                                                                                                |
| T2T CDS genes file                          | Initially empty BED file, populated with chromosome, start, end, and gene_id of the CDS regions in the T2T reference genome during script execution.                                                                                                                  |
| GRCh38 RefSeq GTF file                      | File containing the genes in the GRCh38 reference genome, used in the GRCh38 intersection step.                                                                                                                                                                     |
| GRCh38 CDS genes file                       | Initially empty BED file, populated with chromosome, start, end, and gene_id of the CDS regions in the T2T reference genome during script execution.                                                                                                                 |
| T2T CDS GeneID overlapping only file        | Initially empty BED file, populated with chromosome, start, end, and gene_id of gene_ids that also exist in the GRCh38 CDS file during script execution.                                                                                                             |
| GRCh38 CDS GeneID overlapping only file     | Initially empty BED file, populated with chromosome, start, end, and gene_id of gene_ids that also exist in the T2T CDS file during script execution.                                                                                                                |
| T2T intersected CDS file                    | Initially empty BED file, populated with chromosome, start, end, and gene_id of CDS regions in the T2T CDS GeneID overlapping file that overlap with low coverage CDS regions in the T2T categorical CDS only file during script execution.                          |
| GRCh38 intersected CDS file                 | Initially empty BED file, populated with chromosome, start, end, and gene_id of CDS regions in the GRCh38 CDS GeneID overlapping file that overlap with low coverage CDS regions in the GRCh38 categorical CDS only file during script execution.                      |

Code example:
```
bash low_coverage_CDS_regions.sh T2T_categorical.bed T2T_categorical_CDS_Only.bed GRCh38_categorical.bed GRCh38_categorical_CDS_only.bed T2T_RefSeq_normal_chromosomes.gtf T2T_CDS_genes.bed GRCh38_RefSeq_file.gtf GRCh38_CDS_genes.bed T2T_CDS_GeneID_overlapping.bed GRCh38_CDS_GeneID_overlapping.bed T2T_intersected_CDS.bed GRCh38_intersected_CDS.bed
```

2. Run the [comparing_low_coverage_genes.py](https://github.com/WoutPoelen/Internship_T2T/blob/main/scripts/python/comparing_low_coverage_genes.py) python script like the following code example:
```
python comparing_low_coverage_genes.py T2T_intersected_CDS.bed GRCh38_intersected_CDS.bed shared_genes.txt T2T_unique_genes.txt GRCh38_unqiue_genes.txt
```
3. The plot is saved as Coding_Sequences_in_low_coverage_regions.png
4. The genes which have a coding sequence which overlap with a low coverage region in both files are written to shared_genes.txt. The genes with coding sequences that overlap exclusively with low coverage regions in T2T are recorded in T2T_unique.txt. The genes with coding sequences that overlap exclusively with low coverage regions in GRCh38 are recorded in GRCh38_unique.txt

---
#### Advised next step
Investigate the genes that come out of this step.
