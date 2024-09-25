#!/bin/bash
set -euo pipefail
INFILE=$1
OUTFILE=$2
bcftools query -f "%CHROM %POS %FILTER [%GT] %INFO/SVLEN\n" -i 'INFO/SVTYPE="DEL" || INFO/SVTYPE="INS"' $INFILE > $OUTFILE
