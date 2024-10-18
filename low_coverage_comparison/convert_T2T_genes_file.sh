
# check if a GFF file is provided
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <input_gff_file>"
  exit 1
fi

input_gff=$1
output_gff="converted_${input_gff}"

# Use sed to replace RefSeq chromosome names with standard chromosome names
sed -e 's/NC_060925.1/chr1/g' \
    -e 's/NC_060926.1/chr2/g' \
    -e 's/NC_060927.1/chr3/g' \
    -e 's/NC_060928.1/chr4/g' \
    -e 's/NC_060929.1/chr5/g' \
    -e 's/NC_060930.1/chr6/g' \
    -e 's/NC_060931.1/chr7/g' \
    -e 's/NC_060932.1/chr8/g' \
    -e 's/NC_060933.1/chr9/g' \
    -e 's/NC_060934.1/chr10/g' \
    -e 's/NC_060935.1/chr11/g' \
    -e 's/NC_060936.1/chr12/g' \
    -e 's/NC_060937.1/chr13/g' \
    -e 's/NC_060938.1/chr14/g' \
    -e 's/NC_060939.1/chr15/g' \
    -e 's/NC_060940.1/chr16/g' \
    -e 's/NC_060941.1/chr17/g' \
    -e 's/NC_060942.1/chr18/g' \
    -e 's/NC_060943.1/chr19/g' \
    -e 's/NC_060944.1/chr20/g' \
    -e 's/NC_060945.1/chr21/g' \
    -e 's/NC_060946.1/chr22/g' \
    -e 's/NC_060947.1/chrX/g' \
    -e 's/NC_060948.1/chrY/g' \
    "$input_gff" > "$output_gff"
echo "Chromosome names converted and saved to $output_gff"
