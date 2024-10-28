#!/bin/bash
set -euo pipefail
GENES_INPUT_FILE=$1
SD_INPUT_FILE=$2
CENTROMERES_INPUT_FILE=$3

CDS_OUTPUT_FILE=$4
TRANSCRIPT_OUTPUT_FILE=$5
CENTROMERES_SORTED_OUTPUT_FILE=$6

GRCH38_LOW_COVERAGE_INPUT_FILE=$7
GRCH38_LOW_COVERAGE_CATEGORICAL_FILE=$8

# Checks if there are actually 8 files given. Prints the usage if that is not the case.
if [ "$#" -ne 8 ];
  then
    echo "Usage: genes_GRCh38_file.gtf segmental_duplication.bed centromeres.txt.gz ouput_CDS_bed_file.bed output_transcript.bed output_centromeres.bed T2T_low_coverage_regions.bed T2T_output.bed "
fi

# Processes the Coding sequence file to get the chromosome and start/end location
grep CDS $GENES_INPUT_FILE | cut -f 1,4,5 > $CDS_OUTPUT_FILE
echo $GENES_INPUT_FILE "has been processed and the CDS output is made and send to" $CDS_OUTPUT_FILE

# Processes the transcript file to only get the transcripts and their chromosome, start and end location
awk 'BEGIN{FS="\t"; OFS="\t"} $3 == "transcript" {print $1, $4,$5}' $GENES_INPUT_FILE > $TRANSCRIPT_OUTPUT_FILE
echo $GENES_INPUT_FILE "has been processed and send the transcripts to" $TRANSCRIPT_OUTPUT_FILE

# Processes the CenSat file to get the chromosome, start and end location
cat $CENTROMERES_INPUT_FILE | cut -f 2,3,4 | bedtools sort -i - > $CENTROMERES_SORTED_OUTPUT_FILE
echo $CENTROMERES_INPUT_FILE "has been processed and send the output to" $CENTROMERES_SORTED_OUTPUT_FILE

# Get the amount of times a region overlaps with a category
bedtools intersect -C -a $GRCH38_LOW_COVERAGE_INPUT_FILE -b $CENTROMERES_SORTED_OUTPUT_FILE $TRANSCRIPT_OUTPUT_FILE $CDS_OUTPUT_FILE $SD_INPUT_FILE   -names centromere transcript CDS  SD > $GRCH38_LOW_COVERAGE_CATEGORICAL_FILE
echo $GRCH38_LOW_COVERAGE_INPUT_FILE "has been processed and send the output to" $GRCH38_LOW_COVERAGE_CATEGORICAL_FIL

