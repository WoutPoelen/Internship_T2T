#!/bin/bash
set -euo pipefail
GENES_INPUT_FILE=$1
SD_INPUT_FILE=$2
CENTROMERES_INPUT_FILE=$3

CDS_OUTPUT_FILE=$4
TRANSCRIPT_OUTPUT_FILE=$5
CENTROMERES_SORTED_OUTPUT_FILE=$6

T2T_LOW_COVERAGE_INPUT_FILE=$7
GRCH38_LOW_COVERAGE_INPUT_FILE=$8

T2T_LOW_COVERAGE_CATEGORICAL_FILE=$9
GRCH38_LOW_COVERAGE_CATEGORICAL_FILE=${10}

if [ "$#" -ne 10 ];
  then
    echo "Usage: hg38.ncbiRefSeq.gtf segmental_duplication.bed centromeres.txt.gz ouput_CDS_bed_file
    output_transcript.bed output_centromeres.bed T2T_low_coverage_regions.bed GRCh38_low_coverage regions.bed
    T2T_output.bed GRCh38_output.bed "
fi

# Processes the Coding sequence file to get the chromosome and start/end location
grep CDS $GENES_INPUT_FILE | cut -f 1,4,5 > $CDS_OUTPUT_FILE
echo $GENES_INPUT_FILE " has been processed and the CDS output is made and sent to " $CDS_OUTPUT_FILE

# Processes the transcript file to only get the transcripts and their chromosome, start and end location
awk 'BEGIN{FS="\t"; OFS="\t"} $3 == "transcript" {print $1, $4,$5}' $GENES_INPUT_FILE > $TRANSCRIPT_OUTPUT_FILE
echo $GENES_INPUT_FILE " has been processed and sent the transcripts have been send to " $ TRANSCRIPT_OUTPUT_FILE

# Processes the centromeres file to get the chromosome, start and end location
zcat $CENTROMERES_INPUT_FILE | cut -f 2,3,4 | bedtools sort -i - > $CENTROMERES_SORTED_OUTPUT_FILE
echo "Centromere file has been processed and sent the output to " $CENTROMERES_SORTED_OUTPUT_FILE

# Get the amount of times a region overlaps with a category
bedtools intersect -C -a $GRCH38_LOW_COVERAGE_INPUT_FILE -b $CENTROMERES_SORTED_OUTPUT_FILE $TRANSCRIPT_OUTPUT_FILE $CDS_OUTPUT_FILE $SD_INPUT_FILE   -names centromere transcript CDS  SD > $GRCH38_LOW_COVERAGE_CATEGORICAL_FILE
bedtools intersect -C -a $T2T_LOW_COVERAGE_INPUT_FILE -b $CENTROMERES_SORTED_OUTPUT_FILE $TRANSCRIPT_OUTPUT_FILE $CDS_OUTPUT_FILE $SD_INPUT_FILE   -names centromere transcript CDS  SD > $T2T_LOW_COVERAGE_CATEGORICAL_FILE
echo $GRCH38_LOW_COVERAGE_INPUT_FILE " has been processed and sent the output to " $GRCH38_LOW_COVERAGE_CATEGORICAL_FILE
echo $T2T_LOW_COVERAGE_INPUT_FILE "has been processed and sent the ouput to " $T2T_LOW_COVERAGE_CATEGORICAL_FILE

