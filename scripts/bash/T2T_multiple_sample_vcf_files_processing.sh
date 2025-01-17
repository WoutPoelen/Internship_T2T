#!/bin/bash
non_syntenic_regions=$1
CDS_regions=$2
input_directory=$3
bcftools_output_directory=$4
gene_intersect_output_directory=$5
SV_intersect_output_directory=$6

# Check if directories and files exist
if [ ! -d "$input_directory" ]; then
    echo "Input directory does not exist: $input_directory"
    exit 1
fi

if [ ! -d "$bcftools_output_directory" ]; then
    echo "bcftools output directory does not exist: $bcftools_output_directory"
    exit 1
fi

if [ ! -d "$gene_intersect_output_directory" ]; then
    echo "bedtools output directory does not exist: $gene_intersect_output_directory"
    exit 1
fi

if [ ! -d "$SV_intersect_output_directory" ]; then
    echo "bedtools output directory does not exist: $SV_intersect_output_directory"
    exit 1
fi

if [ ! -f "$non_syntenic_regions" ]; then
    echo "non syntenic file does not exist: $non_syntenic_regions"
    exit 1
fi

if [ ! -f "$CDS_regions" ]; then
    echo "coding sequences file does not exist: $CDS_regions"
    exit 1
fi

# Check if required tools are available
if ! command -v bcftools &> /dev/null; then
    echo "bcftools could not be found"
    exit 1
fi

if ! command -v bedtools &> /dev/null; then
    echo "bedtools could not be found"
    exit 1
fi

# Loops through the samples and obtains the chromosome, start and end position and the type of the structural variation of the structural variations that are larger than 30 base pairs.
# It then writes the structural variations to a directory under the same name as the sample
for file in "$input_directory"/*; do
    # Extract the filename from the path
    filename=$(basename "$file" .vcf).bed
    # Define the output file path
    outfile="$bcftools_output_directory/$filename"
    # Run bcftools query and write to the output file
    bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/SVLEN\t%INFO/SVTYPE\n" -i '(INFO/SVTYPE="DEL" || INFO/SVTYPE="INS") && (INFO/SVLEN > 30 || INFO/SVLEN < -30)' "$file" > "$outfile"
done 

# Intersects CDS regions with the structural variants from the bcftools loop, sorts and remove duplicate CDS, then intersects with non-syntenic regions to identify unique T2T genes overlapping with SVs
for file in "$bcftools_output_directory"/*; do
	filename=$(basename "$file") 
	outfile="$gene_intersect_output_directory/$filename"
	bedtools intersect -a "$CDS_regions" -b "$file" | sort | uniq | bedtools intersect -wb -a - -b "$non_syntenic_regions" > "$outfile"
done

# Intersects CDS regions with non_syntenic and puts into a temporary file to be used later
temp_result=$(mktemp)
bedtools intersect -a "$CDS_regions" -b "$non_syntenic_regions" | sort | uniq > "$temp_result"

# Intersects structural variations from the bcftools loop with the CDS regions in the temporary file to find unique T2T structural variations overlapping with CDS regions
for file in "$bcftools_output_directory"/*; do
	filename=$(basename "$file")
	outfile="$SV_intersect_output_directory/$filename"
	bedtools intersect -a "$file" -b "$temp_result" -wa | sort | uniq > "$outfile"
done

rm "$temp_result"
