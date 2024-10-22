import argparse
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib

def open_file(argument):
    """
    This functions opens the two given bed files that have the coverage of the regions. This file should be made with
    the bash file.

    :param:
        argument.T2T_gz_file (gz file): gz file containing the average coverage of the 500 bp regions from the t2t
        genome.
        argument.GRCh38_gz_file (gz file): gz file containing the average coverage of the 500 bp regions from the GRCh38
        genome.

    :return:
         t2t_dataframe(dataframe): a dataframe that contains the chromosome, start and end of the region.
         Tt also has the coverage of the regions from the t2t bed file.
         hg38_dataframe(dataframe: a dataframe that contains the chromosome, start and end of the region.
         Tt also has the coverage of the regions from the hg38 bed file.
    """
    print("Reading the files")

    # Unzips with gzip and turns the files given by the user into dataframes
    t2t_dataframe = pd.read_csv(argument.T2T_coverage_file, sep="\t", compression="gzip", encoding="utf-8")
    hg38_dataframe = pd.read_csv(argument.GRCh38_coverage_file, sep="\t", compression="gzip", encoding="utf-8")

    # Sets the column names for the dataframe.
    hg38_dataframe.columns = ["Chromosome", "Start", "End", "mean_coverage"]
    t2t_dataframe.columns = ["Chromosome", "Start", "End", "mean_coverage"]

    print("Processing files")

    # Filters the dataframes to get the first 5mbp of the first chromosome. Get rid of this if this isn't necessary
    hg38_dataframe = hg38_dataframe[(hg38_dataframe["End"] <= 5000000) & (hg38_dataframe["Chromosome"] == "chr1")]
    t2t_dataframe = t2t_dataframe[(t2t_dataframe["End"] <= 5000000) & (t2t_dataframe["Chromosome"] == "chr1")]


    return t2t_dataframe, hg38_dataframe

def Counting_the_coverage(dataframe_t2t, dataframe_hg38):
    """
    This functions turns the dataframe column "mean coverage" into a list and makes a dictionary with the amount of
    coverage as the keys and the amount of times that coverage was performed in the file.

    :param:
        dataframe_t2t(dataframe): a dataframe that contains the chromosome, start and end of the region.
        Tt also has the coverage of the regions from the t2t bed file.
        dataframe_hg38(dataframe): a dataframe that contains the chromosome, start and end of the region.
        Tt also has the coverage of the regions from the hg38 bed file.

    :return:
        coverage_occurences_t2t(dictionary): a dictionary that contains the amount of times a amount of base coverage
        occurs in the t2t bed file.
        coverage_occurences_hg38(dictionary): a dictionary that contains the amount of times a amount of base coverage
        occurs in the hg38 bed file.
    """

    # Turns the column "mean_coverage" into a list. This column contains the mean coverage of the bins.
    mean_coverage_t2t = dataframe_t2t["mean_coverage"].tolist()
    mean_coverage_hg38 = dataframe_hg38["mean_coverage"].tolist()

    print("Calculating average coverage occurrences")

    # Calculates amount of times an amount of mean coverage in the files.
    coverage_occurrences_t2t = Counter(mean_coverage_t2t)
    coverage_occurrences_hg38 = Counter(mean_coverage_hg38)

    return coverage_occurrences_t2t, coverage_occurrences_hg38

def calculating_standard_deviation(dataframe_t2t, dataframe_hg38):

    t2t_standard_deviation = dataframe_t2t["mean_coverage"].std()
    hg38_standard_deviation = dataframe_hg38["mean_coverage"].std()


    return t2t_standard_deviation, hg38_standard_deviation


def make_histogram(occurrences_coverage_t2t, occurrences_coverage_hg38):
    """
    Makes the histogram of the lists with the occurrences of the coverage amounts.

    :param:
        occurrences_coverage_t2t(dictionary): a dictionary that contains the amount of times an amount of coverage was
        performed on a base pair.
        occurrences_coverage_hg38 (dictionary): a dictionary that contains the amount of times an amount of coverage was
        performed on a base pair.

    :return:
        saves the plot as coverage_occurrences_histogram.png.
    """
    # Makes the plotting work on the server
    matplotlib.use("Agg")

    fig, ax = plt.subplots()

    print("Generating the plot")

    # Makes a histogram with each file having a different color, automated length of values on the x-axis (bins),
    # a label for the legend, a different alpha value and histtype to make the histograms easier to differentiate
    ax.hist(occurrences_coverage_t2t, color="blue", bins="auto", label="T2T", alpha=0.5, histtype="step")
    ax.hist(occurrences_coverage_hg38,color="green", bins="auto", label="GRCh38", alpha=0.25, histtype="stepfilled")

    # Makes a legend in the upper right
    plt.legend(loc="upper right")

    # Writes the title for the figure in the figure
    plt.title("Mean coverage occurrences")

    # Writes the x and y labels in the figure
    ax.set_xlabel("Mean coverage")
    ax.set_ylabel("Number of occurrences")

    # If a limit for the x and y-axis is necessary comment these two out
    # ax.set_xlim(0, 250)
    # ax.set_ylim(0, 250)

    # Saves the figure in a file
    plt.savefig("coverage_occurrences_histogram.png", bbox_inches="tight")

    print("The average coverage histogram has been generated and is saved as coverage_occurrences_histogram.png.")

    # Makes sure the plot is shown. Commented out because file wouldn't be saved on the server. Get rid of the comment
    # when you want to see the file outside the server
    # plt.show()

def plot_standard_deviation(standard_deviation_t2t, standard_deviation_hg38):
    """
    This function gets the standard deviations of the mean coverage in both reference genomes and plots it.

    :param:
        standard_deviation_t2t (int): standard deviation of the mean coverage column in the T2T dataframe.
        standard_deviation_hg38 (int): standard deviation of the mean coverage column in the GRCh38 dataframe.

    :return:
        The scatterplot is saved as standard_deviation_coverage.png.
    """
    # Makes sure the plot works on the server
    matplotlib.use("Agg")

    print("Generating the plot")

    # Makes a list with values for the y-axis and a list with the reference genomes
    line_values = [standard_deviation_hg38, standard_deviation_t2t]
    line_x_values = ["GRCh38", "T2T"]

    # Sets the figure up with the width and height
    fig, ax = plt.subplots()

    # Creates a scatter plot with the name of the reference genomes in the x-axis and the average mismatch rates in
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

    plt.savefig("standard_deviation_coverage.png", bbox_inches="tight")

    print("The scatterplot has been succesfully generated and saved as standard_deviation_coverage")

    # plt.show()


def main(args):
    t2t_dataframe, hg38_dataframe = open_file(args)
    coverage_occurences_t2t, coverage_occurences_hg38 = Counting_the_coverage(t2t_dataframe, hg38_dataframe)
    t2t_standard_deviation, hg38_standard_deviation = calculating_standard_deviation(t2t_dataframe, hg38_dataframe)
    make_histogram(coverage_occurences_t2t, coverage_occurences_hg38)
    plot_standard_deviation(t2t_standard_deviation, hg38_standard_deviation)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("T2T_coverage_file",
                        help="Path to the file containing the regions of 500 bp and their average coverage "
                             "from the T2T file",
                        metavar="the T2T input gz file")

    parser.add_argument("GRCh38_coverage_file",
                        help="Path to the file containing the regions of 500 bp and their average coverage "
                             "from the GRCh38 file",
                        metavar="the GRCh38 input gz file")

    args = parser.parse_args()
    main(args)
