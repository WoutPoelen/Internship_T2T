#!/bin/bash
set -euo pipefail
INFILE=$1
OUTFILE=$2

# Gets the VCF file from the chosen reference genome and gets the chromosome, position, filter and indels from them
bcftools query -f "%CHROM %POS %FILTER [%GT] %INFO/SVLEN\n" -i 'INFO/SVTYPE="DEL" || INFO/SVTYPE="INS"' $INFILE > $OUTFILE
echo "query complete. The chromosome, position, filter step performed on the sv, whether the sv is heterozygote or
homozygote and the length of the mutation are found in the ${OUTFILE}"
