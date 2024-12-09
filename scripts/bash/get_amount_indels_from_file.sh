#!/bin/bash
set -euo pipefail
INFILE=$1
OUTFILE=$2

# Gets the VCF file from the chosen reference genome and gets the chromosome, position, filter and indels from it
bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/SVLEN\t%INFO/SVTYPE\n" -i '(INFO/SVTYPE="DEL" || INFO/SVTYPE="INS") && (INFO/SVLEN > 30 || INFO/SVLEN < -30)' "$INFILE" > "$OUTFILE"
echo "The chromosome, start position, end position, length and type of the SV has been obtained."
