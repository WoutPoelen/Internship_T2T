#!/bin/bash

non_syntenic_regions="/ifs/data/research/projects/wout/projects/wp1_depth/data/HG002_sample/genes_non_syntenic_regions/chm13v2_non_syntenic_regions.bed"
CDS_regions="/ifs/data/research/projects/wout/projects/wp1_depth/data/HG002_sample/categories_coverage_overlap/comparing_genes/Getting_RefSeq_Only_CDS_Regions/T2T_CDS_RefSeq_only.bed"

input_directory="/ifs/data/research/projects/wout/projects/wp1_depth/data/structural_varation_multiple_samples/GRCh38_VCF_files"
bcftools_output_directory="/ifs/data/research/projects/wout/projects/wp1_depth/data/structural_varation_multiple_samples/bcftools_output_GRCh38"
bedtools_output_directory="/ifs/data/research/projects/wout/projects/wp1_depth/data/structural_varation_multiple_samples/intersect_genes_SV_GRCh38"

# Check if directories exist
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

