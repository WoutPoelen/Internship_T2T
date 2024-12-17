#!/bin/bash
non_syntenic_regions="/ifs/data/research/projects/wout/projects/wp1_depth/data/HG002_sample/genes_non_syntenic_regions/chm13v2_non_syntenic_regions.bed"
CDS_regions="/ifs/data/research/projects/wout/projects/wp1_depth/data/HG002_sample/categories_coverage_overlap/comparing_genes/Getting_RefSeq_Only_CDS_Regions/T2T_CDS_RefSeq_only.bed"

input_directory="/ifs/data/research/projects/wout/projects/wp1_depth/data/structural_varation_multiple_samples/T2T_VCF_files"
bcftools_output_directory="/ifs/data/research/projects/wout/projects/wp1_depth/data/structural_varation_multiple_samples/bcftools_output_T2T"
gene_intersect_output_directory="/ifs/data/research/projects/wout/projects/wp1_depth/data/structural_varation_multiple_samples/non_syntenic_genes"
SV_intersect_output_directory="/ifs/data/research/projects/wout/projects/wp1_depth/data/structural_varation_multiple_samples/non_syntenic_structural_variants"

# Check if directories exist
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

# Loops through the output bed files from the previous for loop and first intersect the coding sequence regions with the structural variant. Then sorts the remaining CDS regions and gets the uniq remaining CDS and intersects the unique CDS regions with the non syntenic regions of T2T to get the unique to T2T genes that overlap with a SV
for file in "$bcftools_output_directory"/*; do
	filename=$(basename "$file") 
	outfile="$gene_intersect_output_directory/$filename"
	bedtools intersect -a "$CDS_regions" -b "$file" | sort | uniq | bedtools intersect -wb -a - -b "$non_syntenic_regions" > "$outfile"
done

# Makes a temporary file to write the output of the next intersect command that intersect the CDS regions and the non syntenic regions, sorts and gets the uniq CDS regions
temp_result=$(mktemp)
bedtools intersect -a "$CDS_regions" -b "$non_syntenic_regions" | sort | uniq > "$temp_result"

# Loops through the output bed files from the first for loop and intersects the structural variations with the CDS regions in the temporary file to get the unique to T2T structural variations that overlap with a CDS region 
for file in "$bcftools_output_directory"/*; do
	filename=$(basename "$file")
	outfile="$SV_intersect_output_directory/$filename"
	bedtools intersect -a "$file" -b "$temp_result" -wa | sort | uniq > "$outfile"
done

rm "$temp_result"
