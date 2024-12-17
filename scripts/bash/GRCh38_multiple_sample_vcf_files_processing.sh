#!/bin/bash

non_syntenic_regions=$1
CDS_regions=$2
input_directory=$3
bcftools_output_directory=$4
bedtools_output_directory=$5

# Check if directories and files exist
if [ ! -d "$input_directory" ]; then
    echo "Input directory does not exist: $input_directory"
    exit 1
fi

if [ ! -d "$bcftools_output_directory" ]; then
    echo "bcftools output directory does not exist: $bcftools_output_directory"
    exit 1
fi

if [ ! -d "$bedtools_output_directory" ]; then
    echo "bedtools output directory does not exist: $bedtools_output_directory"
    exit 1
fi

if [ ! - "$non_syntenic_regions" ]; then
    echo "non syntenic file does not exist: $non_syntenic_regions"
    exit 1
fi

if [ ! -f "$CDS_regions" ]; then
    echo "coding sequence file does not exist: $CDS_regions"
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

# Loops through the output bed files from the previous loop and intersects the coding sequences in GRCh38 with the structural variants obtain in the previous for loop
for file in "$bcftools_output_directory"/*; do
        filename=$(basename "$file")
        outfile="$bedtools_output_directory/$filename"
	bedtools intersect -wb -a "$CDS_regions" -b "$file" | sort | uniq > "$outfile"
done

