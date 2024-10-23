import argparse
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd


def generate_dataframes(argument):
    """
    This functions gets the arguments from the command line, which are a bed file containing average low coverage
    regions from the T2T reference genome, one from the GRCh38 reference genome and one bed file containing the lifted
    over coordinates of the low coverage regions from T2T to GRCh38. After getting the files it turns them into pandas
    dataframes.
    :param:
        argument.T2T_BED_file (bed file): BED file from process_regions_file.py containing the regions with average
        coverage below or equal to 10 from the T2T reference genome.
        argument.GRCh38_BED_file (bed file): BED file from process_regions_file.py containing the regions with average
        coverage below or equal to 10 from the GRCh38 reference genome.
        argument.Liftover_BED_file (bed file): BED file containing the lifted over regions with average coverage below
        or equal to 10 from the T2T reference genome to the GRCh38 reference genome.
    :return:
         t2t_dataframe (dataframe): Pandas dataframe containing the low average coverage BED file regions from the T2T
         reference genome.
         hg38_dataframe (dataframe): Pandas dataframe containing the low average coverage BED file regions from
         the GRCh38 reference genome.
         liftover_dataframe (dataframe): Pandas dataframe containing lifted over low average coverage regions from the
         T2T reference genome to the GRCh38 reference genome.
    """

    # Turns the files gotten from the command line and turns them into pandas dataframes
    t2t_dataframe = pd.read_csv(argument.T2T_BED_file, sep="\t", encoding="utf-8")
    hg38_dataframe = pd.read_csv(argument.GRCh38_BED_file, sep="\t", encoding="utf-8")
    liftover_dataframe = pd.read_csv(argument.Liftover_BED_file, sep="\t", encoding="utf-8")

    print("Reading the files")

    # Adds columns to the dataframes
    t2t_dataframe.columns = ["Chromosome", "Start", "End", "mean_coverage"]
    hg38_dataframe.columns = ["Chromosome", "Start", "End", "mean_coverage"]
    liftover_dataframe.columns = ["Chromosome", "Start", "End", "mean_coverage", "succesful_liftover"]

    print("Processing the files")

    # Filters the pandas frames to only contain the first 5mbp of chromosome 1 and a coverage below 1
    hg38_dataframe = hg38_dataframe[(hg38_dataframe["End"] <= 5000000) & (hg38_dataframe["Chromosome"] == "chr1")]
    t2t_dataframe = t2t_dataframe[(t2t_dataframe["End"] <= 5000000) & (t2t_dataframe["Chromosome"] == "chr1")]
    liftover_dataframe = liftover_dataframe[(liftover_dataframe["End"] <= 5000000) &
                                            (liftover_dataframe["Chromosome"] == "chr1")]

    return t2t_dataframe, hg38_dataframe, liftover_dataframe


def plot_low_coverage(dataframe_t2t, dataframe_hg38, dataframe_liftover):
    """
    This functions plots the dataframes into a scatter plot to visualize overlapping and non-overlapping regions.
    :param:
        t2t_dataframe (dataframe): Pandas dataframe containing the low average coverage BED file regions from the T2T
        reference genome.
        hg38_dataframe (dataframe): Pandas dataframe containing the low average coverage BED file regions from
        the GRCh38 reference genome.
        liftover_dataframe (dataframe): Pandas dataframe containing lifted over low average coverage regions from the
        T2T reference genome to the GRCh38 reference genome.

    :return:
        saves the plot as low_coverage_comparison
    """

    # Makes the plotting work on the server
    matplotlib.use("Agg")

    print("Generating the plot")

    # adds another column to the dataframes showing the origin of the data
    dataframe_t2t["Reference genome"] = "T2T"
    dataframe_hg38["Reference genome"] = "GRCh38"
    dataframe_liftover["Reference genome"] = "liftover_T2T"

    # Plots three scatter plots with the start locations as x values and the origin of the data on the y-axis.
    figs, ax = plt.subplots(figsize=(10, 5))
    plt.scatter(dataframe_t2t["Start"], dataframe_t2t["Reference genome"], s=10)
    plt.scatter(dataframe_liftover["Start"], dataframe_liftover["Reference genome"], s=10)
    plt.scatter(dataframe_hg38["Start"], dataframe_hg38["Reference genome"], s=10)

    # Gives the plot a title
    plt.title("P3-D10 chr1")

    # Gives the x-axis a label
    plt.xlabel("Position (Mbp)")

    # saves the figure into a png file
    plt.savefig("low_coverage_comparison.png")

    print("The low coverage plot has been successfully generated and save as low_coverage_comparison.png")
    # plt.show()


def main(args):
    t2t_dataframe, hg38_dataframe, liftover_dataframe = generate_dataframes(args)
    plot_low_coverage(t2t_dataframe, hg38_dataframe, liftover_dataframe)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("T2T_BED_file",
                        help="Path to the filtered BED file from process_regions_file.py containing "
                             "the regions of 500 bp and the average coverage lower than 10 from the T2T file")
    parser.add_argument("GRCh38_BED_file",
                        help="Path to the filtered BED file from process_regions_file.py containing the regions "
                             "of 500 bp and the average coverage lower than 10 from the GRCh38 file")
    parser.add_argument("Liftover_BED_file",
                        help="Path to the from T2T lifted over BED file containing the coordinates for the locations"
                             "of the lifted over regions on the GRCh38 reference genome")
    args = parser.parse_args()
    main(args)