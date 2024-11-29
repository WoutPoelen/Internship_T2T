#!/bin/bash
set -euo pipefail
INFILE=$1
OUTFILE=$2

# Gets the BAM file from the reference genome and gets the average coverage of regions consisting of 500 base pairs
mosdepth -n --fast-mode --by 500 $OUTFILE $INFILE
echo "Mosdepth analysis complete. Mean coverage of the regions saved to ${OUTFILE}.regions.bed.gz"