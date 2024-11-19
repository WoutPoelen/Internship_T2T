#!/bin/bash
set -euo pipefail

# New input for getting the coding sequences (CDS) regions from the categorical files
T2T_categorical_file=$1
T2T_filtered_categorical_file=$2
GRCh38_categorical_file=$3
GRCh38_filtered_categorical_file=$4

# New input for getting the CDS regions with their gene_id from the RefSeq files
T2T_RefSeq_file=$5
T2T_CDS_genes=$6
GRCh38_RefSeq_file=$7
GRCh38_CDS_genes=$8

# New input for getting the Gene_ids which are in the other CDS genes file
T2T_CDS_GeneID_Overlapping=$9
GRCh38_CDS_GeneID_Overlapping=${10}

# New input for getting the Gene_ids which have CDS that overlap with the low coverage regions:
T2T_intersected_CDS=${11}
GRCh38_intersected_CDS=${12}

# Step 1: Getting the regions that overlap with a CDS region (getting rows with CDS from the categorical file)
awk '$5 == "CDS" && $6 >= 1 {print $1, $2, $3}' OFS='\t' "$T2T_categorical_file" > "$T2T_filtered_categorical_file"
awk '$5 == "CDS" && $6 >= 1 {print $1, $2, $3}' OFS='\t' "$GRCh38_categorical_file" > "$GRCh38_filtered_categorical_file"
echo "The low average coverage regions that overlap with a coding sequence have been obtained and sent to: $T2T_filtered_categorical_file, $GRCh38_filtered_categorical_file"

# Step 2: Getting the CDS regions with their gene_id from the RefSeq files
grep "CDS" "$T2T_RefSeq_file" | awk -F '\t' '{split($9, a, "gene_id \""); split(a[2], b, "\""); print $1, $4, $5, b[1]}' OFS='\t' > "$T2T_CDS_genes"
grep "CDS" "$GRCh38_RefSeq_file" | awk -F '\t' '{split($9, a, "gene_id \""); split(a[2], b, "\""); print $1, $4, $5, b[1]}' OFS='\t' > "$GRCh38_CDS_genes"
echo "The coding sequence rows have been obtained and sent to: $T2T_CDS_genes , $GRCh38_CDS_genes"

# Step 3: Getting the rows with the Gene_id's which are in the other RefSeq file. # The -Fw s done to stop the difference from triggering due to a slight difference.
# Command removes the _1 behind genes on the Y chromosome and of all the other genes. So GRCh38 does it too even though their Y genes don't have a _1.
# This is done to get the genes on the Y chromosomes even though they have a slight name difference. The cut removes the duplicate gene_ID containing the _1
awk '{gsub(/_1$/, "", $4); print $4 "\t" $0}' "$T2T_CDS_genes" | grep -Fwf <(awk '{print $4}' "$GRCh38_CDS_genes"| sort | uniq) | cut -f2- > "$T2T_CDS_GeneID_Overlapping"
awk '{gsub(/_1$/, "", $4); print $4 "\t" $0}' "$GRCh38_CDS_genes" | grep -Fwf <(awk '{print $4}' "$T2T_CDS_genes"| sort | uniq) | cut -f2- > "$GRCh38_CDS_GeneID_Overlapping"
echo "The rows with the Gene_ids which are in both files are gotten and sent to: $T2T_CDS_GeneID_Overlapping , $GRCh38_CDS_GeneID_Overlapping"

# Makes sure that the file is tab seperated
sed -i 's/ \+/\t/g' "$T2T_CDS_GeneID_Overlapping"
sed -i 's/ \+/\t/g' "$GRCh38_CDS_GeneID_Overlapping"

# Step 4: Intersecting the regions from the RefSeq files and their Gene_id's with the low coverage regions found in the
# categorical file.
bedtools intersect -a "$T2T_CDS_GeneID_Overlapping" -b "$T2T_filtered_categorical_file" | bedtools sort -i - | bedtools merge -d 5 -c 4 -o distinct -i - |  uniq > "$T2T_intersected_CDS"
bedtools intersect -a "$GRCh38_CDS_GeneID_Overlapping" -b "$GRCh38_filtered_categorical_file" | bedtools sort -i - | bedtools merge -d 5 -c 4 -o distinct -i - | uniq > "$GRCh38_intersected_CDS"
echo "The regions from the RefSeq files and their Gene_id's with the low coverage regions found in the categorical file are intersected and sent to: $T2T_intersected_CDS , $GRCh38_intersected_CDS"
