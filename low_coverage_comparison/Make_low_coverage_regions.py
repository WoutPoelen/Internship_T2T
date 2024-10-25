import argparse
import pandas as pd

def obtain_files(arguments):
    """
    This functions gets the two gz file from the command line and turns them into dataframes for filtering.

    :param:
        arguments.T2T_gz_file(gz file): gz file containing the 500 bp regions and their coverage mapped to the
        T2T genome.
        arguments.GRCh38_gz_file(gz file): gz file containing the 500 bp regions and their coverage mapped to the
        GRCh38 genome.

    :return:
        t2t_dataframe(dataframe): dataframe containing the 500 bp regions and their coverage from the T2T genome.
        hg38_dataframe(dataframe): dataframe containing the 500 bp regions and their coverage from the GRCh38 genome.
    """
    # Turns the two gz files given in the command line into dataframes and adds column names to them
    t2t_dataframe = pd.read_csv(arguments.T2T_gz_file, sep="\t", compression="gzip", encoding="utf-8")
    hg38_dataframe = pd.read_csv(arguments.GRCh38_gz_file, sep="\t", compression="gzip", encoding="utf-8")
    t2t_dataframe.columns = ["Chromosome", "Start", "End", "mean_coverage"]
    hg38_dataframe.columns = ["Chromosome", "Start", "End", "mean_coverage"]

    return t2t_dataframe, hg38_dataframe


def filtering_function(t2t_dataframe, hg38_dataframe):
    """
    This functions filters the dataframe based on the chromosomes. For each reference genome, one dataframe contains
    rows for the autosomal chromosomes, while another holds rows for the sex chromosomes. This is done, because the
    threshold of the low coverage will be different since the sex chromosomes naturally have fewer reads.

    :param:
        t2t_dataframe(dataframe): dataframe containing the 500 bp regions and their coverage from the T2T genome.
        hg38_dataframe(dataframe): dataframe containing the 500 bp regions and their coverage from the GRCh38 genome.

    :return:
        t2t_dataframe_autosomal(dataframe): dataframe containing the rows from the autosomal chromosomes mapped to
        the T2T genome.
        t2t_dataframe_X_Y(dataframe): dataframe containing the rows from the gendered chromosomes mapped to the T2T
        genome.
        hg38_dataframe_autosomal(dataframe):dataframe containing the rows from the autosomal chromosomes mapped to
        the GRCh38 genome.
        hg38_dataframe_X_Y(dataframe): dataframe containing the rows from the gendered chromosomes mapped to the GRCh38
        genome.
    """

    # List with the chromosomes which will be used for the median calculations
    autosomal_chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                             "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                             "chr20", "chr21", "chr22"]

    # List with the gendered chromosomes which will not be used for the median calculations, but will be used for
    # comparisons between the two reference genomes
    genedered_chromosomes = ["chrX", "chrY"]

    # Putting the chromosomes and their respective rows in the correct dataframe. Technique of filtering are different,
    # because GRCh38 contains "alternative contigs" which cannot be removed with the T2T lines (chr1_KI270713v1_random,
    # chrY_KI270740v1_random)
    t2t_dataframe_autosomal = t2t_dataframe[t2t_dataframe["Chromosome"].str.contains("chrX|chrY|chrM") == False]
    t2t_dataframe_X_Y = t2t_dataframe[t2t_dataframe["Chromosome"].str.contains("chrX|chrY")]

    hg38_dataframe_autosomal = hg38_dataframe[hg38_dataframe["Chromosome"].isin(autosomal_chromosomes)]
    hg38_dataframe_X_Y = hg38_dataframe[hg38_dataframe["Chromosome"].isin(genedered_chromosomes)]


    return t2t_dataframe_autosomal, t2t_dataframe_X_Y ,hg38_dataframe_autosomal, hg38_dataframe_X_Y


def getting_low_coverage_regions(t2t_dataframe_autosomal, t2t_dataframe_X_Y,
                                 hg38_dataframe_autosomal, hg38_dataframe_X_Y):
    """
    This function gets the rows from the autosomal chromosomes dataframes which have an average coverage lower than
    1/3rd of the median from their respective reference genome. It also gets the rows from the gendered chromosomes
    dataframes which have an average coverage lower than 1/6th of the median.

    :param:
        t2t_dataframe_autosomal(dataframe): dataframe containing the rows from the autosomal chromosomes mapped to
        the T2T genome.
        t2t_dataframe_X_Y(dataframe): dataframe containing the rows from the gendered chromosomes mapped to the T2T
        genome.
        hg38_dataframe_autosomal(dataframe): hg38_dataframe_autosomal(dataframe):dataframe containing the rows from the
        autosomal chromosomes mapped to the GRCh38 genome.
        hg38_dataframe_X_Y(dataframe): dataframe containing the rows from the gendered chromosomes mapped to the GRCh38
        genome.

    :return:
        autosomal_t2t_dataframe_filtered(dataframe): dataframe containing the rows from the autosomal chromosomes with
        a coverage lower than 1/3rd of the median.
        autosomal_hg38_dataframe_filtered(dataframe): dataframe containing the rows from the autosomal chromosomes with
        a coverage lower than 1/3rd of the median.
        x_y_t2t_dataframe_filtered(dataframe): dataframe containing the rows from the gendered chromosomes with
        a coverage lower than 1/6th of the median.
        x_y_hg38_dataframe_filtered(dataframe): dataframe containing the rows from the gendered chromosomes with
        a coverage lower than 1/6th of the median.
    """
    # Calculates the median of the average coverages of the 500 bp regions in autosomal chromosomes and dividing it
    # by three to get the threshold of what is considered low coverage for these chromosomes
    threshold_autosomal_t2t = t2t_dataframe_autosomal["mean_coverage"].median() / 1/3
    threshold_autosomal_hg38 = hg38_dataframe_autosomal["mean_coverage"].median() / 1/3

    # Calculates the median of the average coverages of the 500 bp regions in the X and Y chromosomes and dividing it
    # by six to get the threshold of what is conidered low coverage for these chromosomes
    threshold_X_Y_t2t = t2t_dataframe_X_Y["mean_coverage"].median() / 1/6
    threshold_X_y_hg38 = hg38_dataframe_X_Y["mean_coverage"].median() / 1/6

    # Different calculations were used, because the X and Y chromosomes naturally get fewer reads due to HG002
    # having one X and one Y chromosome

    # Prints the thresholds of the low coverage
    print("The low coverage threshold for the autosomal chromosomes in T2T is:", threshold_autosomal_t2t, "\n",
          "The low coverage threshold for the autosomal chromosomes in GRCh38 is: ", threshold_autosomal_hg38, "\n",
          "The low coverage threshold for the X and Y chromosomes in T2T is:", threshold_X_Y_t2t, "\n",
          "The low coverage threshold for the X and Y chromosomes in GRCh38 is:", threshold_X_y_hg38)

    # Gets the rows of each dataframe which have a lower average coverage than their respective threshold
    autosomal_t2t_dataframe_filtered = t2t_dataframe_autosomal[t2t_dataframe_autosomal["mean_coverage"]
                                                               < threshold_autosomal_t2t]
    x_y_t2t_dataframe_filtered = t2t_dataframe_X_Y[t2t_dataframe_X_Y["mean_coverage"] < threshold_X_Y_t2t]

    autosomal_hg38_dataframe_filtered = hg38_dataframe_autosomal[hg38_dataframe_autosomal["mean_coverage"]
                                                                < threshold_autosomal_hg38]
    x_y_hg38_dataframe_filtered = hg38_dataframe_X_Y[hg38_dataframe_X_Y["mean_coverage"] < threshold_X_y_hg38]


    return (autosomal_t2t_dataframe_filtered, autosomal_hg38_dataframe_filtered, x_y_t2t_dataframe_filtered,
            x_y_hg38_dataframe_filtered)

def write_to_files(autosomal_t2t_dataframe_filtered, autosomal_hg38_dataframe_filtered, x_y_t2t_dataframe_filtered,
                   x_y_hg38_dataframe_filtered, arguments):
    """
    This function combines the autosomal and X/Y dataframes for each reference genome and writes it to the files (one
    for T2T and one for GRCh38) which are given in the command line.

    :param:
        t2t_dataframe_autosomal(dataframe): dataframe containing the rows from the autosomal chromosomes mapped to
        the T2T genome.
        t2t_dataframe_X_Y(dataframe): dataframe containing the rows from the gendered chromosomes mapped to the T2T
        genome.
        hg38_dataframe_autosomal(dataframe): hg38_dataframe_autosomal(dataframe):dataframe containing the rows from the
        autosomal chromosomes mapped to the GRCh38 genome.
        hg38_dataframe_X_Y(dataframe): dataframe containing the rows from the gendered chromosomes mapped to the GRCh38
        genome.


    :return:
        arguments.low_coverage_output_T2T(bed file): an empty bed file where the combined t2t dataframe will be
        written to.
        arguments.low_coverage_output_GRCh38(bed file): an empty bed file where the combined GRCh38 dataframe will be
        written to.
    """
    # Combines the T2T and GRCh38 dataframes. So they can be written to empty bed files
    combined_t2t_dataframe = pd.concat([autosomal_t2t_dataframe_filtered, x_y_t2t_dataframe_filtered],
                                       ignore_index=True)
    combined_hg38_dataframe = pd.concat([autosomal_hg38_dataframe_filtered, x_y_hg38_dataframe_filtered],
                                        ignore_index=True)

    # Writes the dataframes to their given bed files
    combined_t2t_dataframe.to_csv(arguments.low_coverage_output_T2T, sep="\t", index=False, header=False)
    combined_hg38_dataframe.to_csv(arguments.low_coverage_output_GRCh38, sep="\t", index=False, header=False)

    print("The T2T dataframe with the low coverage regions has been written to:", arguments.low_coverage_output_T2T,
          "\n", "The GRCh38 dataframe with the low coverage regions has been written to:",
          arguments.low_coverage_output_GRCh38 )

def main(args):
    t2t_dataframe, hg38_dataframe = obtain_files(args)

    t2t_dataframe_autosomal, t2t_dataframe_X_Y ,hg38_dataframe_autosomal, hg38_dataframe_X_Y = (
        filtering_function(t2t_dataframe, hg38_dataframe))

    (autosomal_t2t_dataframe_filtered, autosomal_hg38_dataframe_filtered, x_y_t2t_dataframe_filtered,
     x_y_hg38_dataframe_filtered) = getting_low_coverage_regions(t2t_dataframe_autosomal, t2t_dataframe_X_Y,
                                                                 hg38_dataframe_autosomal, hg38_dataframe_X_Y)

    write_to_files(autosomal_t2t_dataframe_filtered, autosomal_hg38_dataframe_filtered, x_y_t2t_dataframe_filtered,
            x_y_hg38_dataframe_filtered, args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("T2T_gz_file",
                        help="Path to the gz file containing the 500 bp regions and their average coverage of the T2T"
                             "reference genome.")
    parser.add_argument("GRCh38_gz_file",
                        help="Path to the gz file containing the 500 bp regions and their average coverage of the"
                             "GRCh38 reference genome.")
    parser.add_argument("low_coverage_output_T2T",
                        help="an empty bed file where the filtered t2t dataframes will be written to")
    parser.add_argument("low_coverage_output_GRCh38",
                        help="an empty bed file where the filtered GRCh38 dataframes will be written to")
    args = parser.parse_args()
    main(args)