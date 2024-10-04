import argparse
import pandas as pd

def openfile(argument):
    """
    Gets the input files as arguments from the server and turns them into dataframes
    :param:
        arg T2T_BED_file (gz file): a gz BED file containing the regions of the T2T genome that needs to be uncompressed
        arg GRCh38_BED_file (gz file): a gz BED file containing the regions of the GRCh38 genome which needs
        to be uncompressed

    :return:
        t2t_dataframe (dataframe): a dataframe containing the uncompressed regions of the T2T genome.
        hg38_dataframe (dataframe): a dataframe containing the uncompressed regions of the GRCh38 genome.
    """

    # Gets the designated gz file and decompresses it with gzip, decodes it with utf-8 and makes it into a dataframe
    t2t_dataframe = pd.read_csv(argument.T2T_GZ_file, sep="\t", compression="gzip", encoding="utf-8")
    hg38_dataframe = pd.read_csv(argument.GRCh38_GZ_file, sep="\t", compression="gzip", encoding="utf-8")


    return t2t_dataframe, hg38_dataframe

def process_dataframes(dataframe_t2t, dataframe_hg38):
    """
    Gives column names to the dataframes and filters the dataframe on the subset of 5mbp, chromosome 1 and the
    average coverage which needs to be below or exactly 10.
    :param:
        dataframe_t2t (dataframe): a dataframe containing the uncompressed regions of the T2T genome.
        dataframe_hg38 (dataframe: a dataframe containing the uncompressed regions of the GRCh38 genome.
    :return:
         t2t_dataframe (dataframe): a dataframe containing the regions on the 5mbp chromosome 1 subset with a mean
         coverage below or exactly 10.
         hg38_dataframe (dataframe): a dataframe containing the regions on the 5bmp chromosome 1 subset with a mean
         coverage below or exactly 10.
    """

    # Adds the column names to the dataframes
    dataframe_t2t.columns = ["Chromosome", "Start", "End", "mean_coverage"]
    dataframe_hg38.columns = ["Chromosome", "Start", "End", "mean_coverage"]

    # Filters the dataframe on the location of the regions, which needs below the 5mbp and on chromosome 1 and the mean
    # coverage of the region needs to be below 10. Remove the location filter if the entire genome needs to be done.
    hg38_dataframe = dataframe_hg38[(dataframe_hg38["End"] <= 5000000) & (dataframe_hg38["Chromosome"] == "chr1") &
                                    (dataframe_hg38["mean_coverage"] <= 10)]
    t2t_dataframe = dataframe_t2t[(dataframe_t2t["End"] <= 5000000) & (dataframe_t2t["Chromosome"] == "chr1") &
                                  (dataframe_t2t["mean_coverage"] <= 10)]

    return t2t_dataframe, hg38_dataframe

def write_to_file(t2t_dataframe, hg38_dataframe, argument):
    """

    :param:
        t2t_dataframe (dataframe): a dataframe containing the regions on the 5mbp chromosome 1 subset with a mean
        coverage below or exactly 10.
        hg38_dataframe (dataframe): a dataframe containing the regions on the 5bmp chromosome 1 subset with a mean
        coverage below or exactly 10.
        arg output_file:
    :return:
        arg.output_file_T2T (BED file): a BED file containing the filtered regions from the t2t regions.
        arg.output_file_GRCh38 (BED file): a BED file containing the regions of the GRCh38 genome.
    """

    # Writes the dataframe without the headers and indexes to the argument output file.
    # Newline causes the extra newlines are not written.
    with open(argument.output_file_T2T, "w", newline="") as T2T_file:
        t2t_dataframe.to_csv(T2T_file, sep="\t", index=False, header=False)

    # Writes the dataframe without the headers and indexes to the argument output file.
    # Newline causes the extra newlines are not written.
    with open(argument.output_file_GRCh38, "w", newline="") as HG38_file:
        hg38_dataframe.to_csv(HG38_file, sep="\t", index=False, header=False)


def main(args):
    t2t_dataframe, hg38_dataframe = openfile(args)
    t2t_dataframe, hg38_dataframe = process_dataframes(t2t_dataframe, hg38_dataframe)
    write_to_file(t2t_dataframe, hg38_dataframe, args)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Turns two gz file (one for T2T and one for GRCh38) into a BED file ")
    parser.add_argument("T2T_GZ_file",
                        help="Path to GZ BED file containing the regions of 500 bp and the average coverage from the T2T file",
                        metavar="T2T BED GZ input file")
    parser.add_argument("GRCh38_GZ_file",
                        help="Path to GZ BED file containing the regions of 500 bp and the average coverage from the GRCh38 file",
                        metavar="GRCh38 GZ BED input file")
    parser.add_argument("output_file_T2T",
                        help="Path to BED output file in which the T2T dataframe is written",
                        metavar="T2T BED output file")
    parser.add_argument("output_file_GRCh38",
                        help="Path to BED output file in which the GRCh38 dataframe is written",
                        metavar="GRCh38 BED output file")
    args = parser.parse_args()
    main(args)


# Alternative aproach:

#!/bin/bash
# set -euo pipefail
IN_GZ=$1
OUT_BED=$2
# zcat ${IN_GZ} | awk '{FS="\t"; OFS="\t"} (($3 < 500000) && ($1 == 'chr1')) {print $1,$2,$3,$5(?or which ever has the depth)}' > ${OUT_BED}

