import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib


def open_file(argument):
    """
    This functions opens the two given bed files that have the coverage of the regions. These files should be made with
    the bash file get_low_coverage_regions.sh.

    :param:
        argument.T2T_coverage_file (gz file): gz file containing the average coverage of the 500 bp regions from the t2t
        file.
        argument.GRCh38_coverage_file (gz file): gz file containing the average coverage of the 500 bp regions from
        the GRCh38 file.

    :return:
         t2t_dataframe(dataframe): a dataframe that contains the chromosome, start and end of the region.
         It also has the coverage of the regions from the t2t bed file.
         hg38_dataframe(dataframe): a dataframe that contains the chromosome, start and end of the region.
         It also has the coverage of the regions from the hg38 bed file.
    """
    print("Reading the files")
    column_names = ["Chromosome", "Start", "End", "mean_coverage"]

    # Unzips with gzip and turns the files given by the user into dataframes
    t2t_dataframe = pd.read_csv(argument.T2T_coverage_file, sep="\t", compression="gzip", encoding="utf-8",
                                names=column_names)
    hg38_dataframe = pd.read_csv(argument.GRCh38_coverage_file, sep="\t", compression="gzip", encoding="utf-8",
                                 names=column_names)

    autosomal_chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                             "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                             "chr20", "chr21", "chr22", "chrX", "chrY"]

    hg38_dataframe = hg38_dataframe[hg38_dataframe["Chromosome"].isin(autosomal_chromosomes)]
    t2t_dataframe = t2t_dataframe[t2t_dataframe["Chromosome"].isin(autosomal_chromosomes)]

    return t2t_dataframe, hg38_dataframe

def get_statistics(t2t_dataframe, hg38_dataframe):
    """
    This function calculates the standard deviation and mean of the coverage of the regions from the t2t and hg38
    bed file. It also calculates the amount of regions with a coverage above a certain threshold.

    :param:
        t2t_dataframe(dataframe): a dataframe that contains the chromosome, start and end of the region.
        It also has average coverage of the regions from the t2t bed file.
        hg38_dataframe(dataframe): a dataframe that contains the chromosome, start and end of the region.
        It also has the average coverage of the regions from the hg38 bed file.

    :return:
        t2t_standard_deviation(float): the standard deviation of the coverage from the t2t bed file.
        hg38_standard_deviation(float): the standard deviation of the coverage from the hg38 bed file.
        t2t_mean(float): the mean of the coverage from the t2t file.
        hg38_mean(float): the mean of the coverage from the hg38 file.
    """
    # Calculates the amount of regions that have a coverage higher than 80. Change the number depending on what the
    # threshold is meant to be
    t2t_high_coverage = len(t2t_dataframe[t2t_dataframe["mean_coverage"] >= 80])
    hg38_high_coverage = len(hg38_dataframe[hg38_dataframe["mean_coverage"] >= 80])

    # Calculates the standard deviation of the mean coverage
    t2t_standard_deviation = t2t_dataframe["mean_coverage"].std()
    hg38_standard_deviation = hg38_dataframe["mean_coverage"].std()
    print("T2T standard deviation:", t2t_standard_deviation)
    print("HG38 standard deviation:", hg38_standard_deviation)

    # Calculates the mean of the mean coverage
    t2t_mean = t2t_dataframe["mean_coverage"].mean()
    hg38_mean = hg38_dataframe["mean_coverage"].mean()
    print("T2T mean:", t2t_mean)
    print("HG38 mean:", hg38_mean)

    return t2t_standard_deviation, hg38_standard_deviation, t2t_mean, hg38_mean, t2t_high_coverage, hg38_high_coverage


def make_histogram( t2t_standard_deviation, hg38_standard_deviation, t2t_mean, hg38_mean,
                   t2t_dataframe, hg38_dataframe, t2t_high_coverage, hg38_high_coverage):
    """
    Makes the histogram of the mean coverage columns. Also plots the mean and standard deviation of the coverage.

    :param:
        t2t_standard_deviation(float): the standard deviation of the coverage from the t2t bed file.
        hg38_standard_deviation(float): the standard deviation of the coverage from the hg38 bed file.
        t2t_mean(float): the mean of the coverage from the t2t file.
        hg38_mean(float): the mean of the coverage from the hg38 file.
        t2t_dataframe(dataframe): a dataframe that contains the chromosome, start and end of the region.
        It also has the average coverage of the regions in the t2t bed file.
        hg38_dataframe(dataframe): a dataframe that contains the chromosome, start and end of the region.
        It also has the average coverage of the regions from the hg38 bed file.

    :return:
        saves the plot as coverage_occurrences_histogram.png.
    """
    # Makes the plotting work on the server
    matplotlib.use("Agg")

    fig, ax = plt.subplots()

    # Changes the average coverage columns of the dataframes into lists
    data_t2t = t2t_dataframe["mean_coverage"].tolist()
    data_hg38 = hg38_dataframe["mean_coverage"].tolist()

    print("Generating the coverage occurrence histogram")

    # Makes the bin_size the exactly one by going to the highest number out of both lists
    bin_size = np.arange(0, max(max(data_t2t), max(data_hg38)))

    # Makes a histogram with each file having a different color, a bin size of one,
    # a label for the legend, a different alpha value and histtype to make the histograms easier to differentiate
    ax.hist(data_t2t, color="blue", bins=bin_size, label="T2T", alpha=0.5, histtype="step")
    ax.hist(data_hg38, color="green", bins=bin_size, label="GRCh38", alpha=0.25, histtype="stepfilled")

    # Adds the high coverage regions to the histogram
    plt.bar(80, height=t2t_high_coverage, width=1, color="blue", alpha=0.5)
    plt.bar(80, height=hg38_high_coverage, width=1, color="green", alpha=0.25)

    # Makes a legend in the upper right
    plt.legend(loc="upper right", prop={"size": 8})

    # Writes the title for the figure in the figure
    plt.title("Mean coverage occurrences")

    # Writes the x and y labels in the figure
    ax.set_xlabel("Mean coverage")
    ax.set_ylabel("Number of occurrences")

    # The x-axis isn't going past the 100. Otherwise, the plot would have unnecessary numbers
    ax.set_xlim(0, 80)

    # Changes the 100 x tick to 100+ to show that there are still regions with coverages higher than 10
    xticks = ax.get_xticks()
    xtick_labels = [str(int(tick))
                    if tick != 80 else "80+"
                    for tick in xticks]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels)

    # Saves the figure in a file
    plt.savefig("coverage_occurrences_histogram.png", bbox_inches="tight")

    print("The average coverage histogram has been generated and is saved as coverage_occurrences_histogram.png.")

    # Makes sure the plot is shown. Commented out because file wouldn't be saved on the server. Get rid of the comment
    # when you want to see the file outside the server
    # plt.show()

def plot_standard_deviation(standard_deviation_t2t, standard_deviation_hg38):
    """
    This function gets the standard deviations of the mean coverage in both reference genomes and plots it in a
    scatterplot.

    :param:
        standard_deviation_t2t (float): standard deviation of the mean coverage column in the T2T dataframe.
        standard_deviation_hg38 (float): standard deviation of the mean coverage column in the GRCh38 dataframe.

    :return:
        The scatterplot is saved as standard_deviation_coverage.png.
    """
    # Makes sure the plot works on the server
    matplotlib.use("Agg")

    print("Generating the standard deviation scatterplot")

    # Makes a list with values for the y-axis and a list with the reference genomes
    line_values = [standard_deviation_hg38, standard_deviation_t2t]
    line_x_values = ["GRCh38", "T2T"]

    # Sets the figure up with the width and height
    fig, ax = plt.subplots()

    # Creates a scatter plot with the name of the reference genomes on the x-axis and the standard deviation on
    # the y-axis in the figure
    ax.scatter(x="GRCh38",y= standard_deviation_hg38, color="blue", label="GRCh38")
    ax.scatter(x = "T2T", y = standard_deviation_t2t, color="green", label="T2T" )
    ax.plot(line_x_values, line_values, color="red", linestyle="dashed", label="Difference")

    # Sets a limit to the x-axis and y-axis and the label for the y-axis
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylabel("Coverage standard deviation")

    # Sets the title and legend up
    plt.title("Standard deviation of the mean coverage")
    plt.legend()

    # Saves the plot as standard deviation_coverage.png
    plt.savefig("standard_deviation_coverage.png", bbox_inches="tight")

    print("The scatterplot has been succesfully generated and saved as standard_deviation_coverage")

    # plt.show()


def main(args):
    t2t_dataframe, hg38_dataframe= open_file(args)

    t2t_standard_deviation, hg38_standard_deviation, t2t_mean, hg38_mean, t2t_high_coverage, hg38_high_coverage\
        = get_statistics(t2t_dataframe, hg38_dataframe)

    make_histogram(t2t_standard_deviation, hg38_standard_deviation, t2t_mean, hg38_mean, t2t_dataframe, hg38_dataframe,
                   t2t_high_coverage, hg38_high_coverage)

    plot_standard_deviation(t2t_standard_deviation, hg38_standard_deviation)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("T2T_coverage_file",
                        help="Path to the gz file containing the regions of 500 bp and their average coverage "
                             "from the T2T file",
                        metavar="the T2T input gz file")

    parser.add_argument("GRCh38_coverage_file",
                        help="Path to the gz file containing the regions of 500 bp and their average coverage "
                             "from the GRCh38 file",
                        metavar="the GRCh38 input gz file")

    args = parser.parse_args()
    main(args)
