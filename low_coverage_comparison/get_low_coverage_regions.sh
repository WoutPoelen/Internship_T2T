#!/bin/bash
set -euo pipefail

T2T_GENOME_REGIONS_INPUT_FILE=$1
T2T_LOW_COVERAGE_REGIONS_OUTPUT=$2

GRCh38_GENOME_REGIONS_INPUT=$3
GRCh38_LOW_COVERAGE_REGION_OUTPUT=$4

# Processes the from Mosdepth obtained T2T gz file and gets the regions which have a average coverage lower than 10
zcat $T2T_GENOME_REGIONS_INPUT_FILE | awk '{FS="\t";OFS="\t"} ($4 <10) {print $1,$2,$3,$4}' > $T2T_LOW_COVERAGE_REGIONS_OUTPUT
echo $T2T_GENOME_REGIONS_INPUT_FILE "has been processed and the output has been send to" $T2T_LOW_COVERAGE_REGIONS_OUTPUT

# Processes the from Mosdepth obtained GRCh38 gz file and gets the regions which have a average coverage lower than 10
zcat $GRCh38_GENOME_REGIONS_INPUT | awk '{FS="\t";OFS="\t"} ($4 <10) {print $1,$2,$3,$4}' > $GRCh38_LOW_COVERAGE_REGION_OUTPUT
echo GRCh38_GENOME_REGIONS_INPUT "has been processed and the output has been send to" $GRCh38_LOW_COVERAGE_REGION_OUTPUT

