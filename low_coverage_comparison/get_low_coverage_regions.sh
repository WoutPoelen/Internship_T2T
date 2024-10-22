#!/bin/bash
set -euo pipefail

T2T_GENOME_REGIONS_INPUT_FILE=$1
T2T_LOW_COVERAGE_REGIONS_OUTPUT=$2

GRCh38_GENOME_REGIONS_INPUT=$3
GRCh38_LOW_COVERAGE_REGION_OUTPUT=$4

# Calculates the median
GRCh38_median=$(zcat "$GRCh38_GENOME_REGIONS_INPUT" | sort -n | awk '{arr[NR]=$4} END {if (NR%2) print arr[(NR+1)/2]; else print (arr[NR/2]+arr[NR/2+1])/2}')
echo "The median of the GRCh38 file is: $GRCh38_median"
T2T_median=$(zcat "$T2T_GENOME_REGIONS_INPUT_FILE" | sort -n | awk '{arr[NR]=$4} END {if (NR%2) print arr[(NR+1)/2]; else print (arr[NR/2]+arr[NR/2+1])/2}')
echo "The median of the T2T file is: $T2T_median"

# Processes the from Mosdepth obtained T2T gz file and gets the regions which have a average coverage lower than the median divided by 3
zcat "$T2T_GENOME_REGIONS_INPUT_FILE" | awk -v median_t2t="$T2T_median" '{FS="\t";OFS="\t"} ($4 < (median_t2t / 3)) {print $1,$2,$3,$4}' > "$T2T_LOW_COVERAGE_REGIONS_OUTPUT"
echo "$T2T_GENOME_REGIONS_INPUT_FILE has been processed and the output has been send to $T2T_LOW_COVERAGE_REGIONS_OUTPUT"

# Processes the from Mosdepth obtained GRCh38 gz file and gets the regions which have a average coverage lower than the median divided by 3
zcat "$GRCh38_GENOME_REGIONS_INPUT" | awk -v median_hg38="$GRCh38_median"'{FS="\t";OFS="\t"} ($4 < (median_hg38 / 3)) {print $1,$2,$3,$4}' > "$GRCh38_LOW_COVERAGE_REGION_OUTPUT"
echo "$GRCh38_GENOME_REGIONS_INPUT has been processed and the output has been send to $GRCh38_LOW_COVERAGE_REGION_OUTPUT"

