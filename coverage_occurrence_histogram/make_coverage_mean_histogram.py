import gzip
import os
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt

def open_file():
    """
    This functions opens the two given bed files that have the coverage of the regions. This file should be made with
    the bash file:

    Args:
    nothing

    Returns:
         t2t_dataframe(dataframe): a dataframe that contains the chromosome, start and end of the region.
         Tt also has the coverage of the regions from the t2t bed file.
         hg38_dataframe(dataframe: a dataframe that contains the chromosome, start and end of the region.
         Tt also has the coverage of the regions from the hg38 bed file.
    """
    # Asks for the input path for the t2t bed regions file
    input_path_t2t = input("Give the path to the t2t bed regions file: ")
    print("given path: " + input_path_t2t)

    # Checks if there is a file in the file path
    if os.path.exists(input_path_t2t):
        print("File exists")
        print(" ")
        # Opens the file with gunzip and turns it into a pandas dataframe
        with gzip.open(input_path_t2t, "r") as f:
            t2t_dataframe = pd.read_csv(input_path_t2t, sep="\t")

    else:
        print("Nothing found in the file path")

    # Asks for the input path for the t2t BAM file
    input_path_hg38 = input("Give the path to the GRCh38 bed regions file: ")
    print("given path: " + input_path_hg38)

    # Checks if there is a file in the file path
    if os.path.exists(input_path_hg38):
        print("File exists")
        print(" ")
        # Opens the file with gunzip and turns it into a dataframe
        with gzip.open(input_path_hg38, "r") as f:
            hg38_dataframe = pd.read_csv(input_path_hg38, sep="\t")


    else:
        print("Nothing found in the file path")

    # Sets the column names for the dataframe.
    hg38_dataframe.columns = ["Chromosome", "Start", "End", "mean_coverage"]
    t2t_dataframe.columns = ["Chromosome", "Start", "End", "mean_coverage"]


    return t2t_dataframe, hg38_dataframe

def Counting_the_coverage(dataframe_t2t, dataframe_hg38):
    """
    This functions turns the dataframe column "mean coverage" into a list and makes a dictionary with the amount of
    coverage as the keys and the amount of times that coverage was performed in the file.

    Args:
        dataframe_t2t(dataframe): a dataframe that contains the chromosome, start and end of the region.
        Tt also has the coverage of the regions from the t2t bed file.
        dataframe_hg38(dataframe): a dataframe that contains the chromosome, start and end of the region.
        Tt also has the coverage of the regions from the hg38 bed file.

    Returns:
        coverage_occurences_t2t(dictionary): a dictionary that contains the amount of times a amount of base coverage
        occurs in the t2t bed file.
        coverage_occurences_hg38(dictionary): a dictionary that contains the amount of times a amount of base coverage
        occurs in the hg38 bed file.
    """
    # Turns the column "mean_coverage" into a list. This column contains the mean coverage of the bins.
    mean_coverage_t2t = dataframe_t2t["mean_coverage"].tolist()
    mean_coverage_hg38 = dataframe_hg38["mean_coverage"].tolist()

    # Calculates amount of times an amount of mean coverage in the files.
    coverage_occurences_t2t = Counter(mean_coverage_t2t)
    coverage_occurences_hg38 = Counter(mean_coverage_hg38)

    return coverage_occurences_t2t, coverage_occurences_hg38


def make_histogram(occurrences_coverage_t2t, occurrences_coverage_hg38):
    """
    Makes the histogram of the lists with the occurrences of the coverage amounts.

    Args:
        occurrences_coverage_t2t(dictionary): a dictionary that contains the amount of times an amount of coverage was
        performed on a base pair.
        occurrences_coverage_hg38 (dictionary): a dictionary that contains the amount of times an amount of coverage was
        performed on a base pair.

    :return:
        saves the plot as coverage_occurrences_histogram.jpg.
    """
    # Makes a histogram with each file having a different color, automated length of values on the x-axis (bins),
    # a label for the legend, a different alpha value and histtype to make the histograms easier to differentiate
    plt.hist(occurrences_coverage_t2t, color="blue", bins="auto", label="T2T", alpha=0.5, histtype="step")
    plt.hist(occurrences_coverage_hg38,color="green", bins="auto", label="GRCh38", alpha=0.25, histtype="stepfilled")

    # Makes a legend in the upper right
    plt.legend(loc="upper right")

    # Writes the title for the figure in the figure
    plt.title("Mean coverage occurences")

    # Writes the x and y labels in the figure
    plt.xlabel("Mean coverage")
    plt.ylabel("Number of occurrences")

    # Saves the figure in a file
    plt.savefig("coverage_occurrences_histogram.jpg", bbox_inches="tight")

    # Makes sure the plot is shown
    plt.show()


if __name__ == "__main__":
    t2t_dataframe, hg38_dataframe = open_file()
    coverage_occurences_t2t, coverage_occurences_hg38 = Counting_the_coverage(t2t_dataframe, hg38_dataframe)
    make_histogram(coverage_occurences_t2t, coverage_occurences_hg38)
