#!/bin/bash
set -euo pipefail
INFILE=$1
OUTFILE=$2
mosdepth -n --fast-mode --by 500 $OUTFILE $INFILE
echo "Mosdepth analysis complete. Mean coverage of the regions saved to ${OUTFILE}.regions.bed.gz"