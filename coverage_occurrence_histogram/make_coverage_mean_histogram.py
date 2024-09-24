import gzip
import os
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt

# T2T: C:\Users\Z659216\Desktop\BAM_files\test_run2_T2T.regions.bed.gz
# Hg38: C:\Users\Z659216\Desktop\BAM_files\test_run_Hg38.regions.bed.gz

def open_file():


    input_path_t2t = input("Give the path to the t2t bed regions file: ")
    print("given path: " + input_path_t2t)

    # Checks if there is a file in the file path
    if os.path.exists(input_path_t2t):
        print("File exists")
        print(" ")
    #     Opens the file and turns it into a pandas dataframe.
        with gzip.open(input_path_t2t, "r") as f:
            t2t_dataframe = pd.read_csv(input_path_t2t, sep="\t")

    else:
        print("Nothing found in the file path")


    input_path_hg38 = input("Give the path to the GRCh38 bed regions file: ")
    print("given path: " + input_path_hg38)

    # Checks if there is a file in the file path
    if os.path.exists(input_path_hg38):
        print("File exists")
        print(" ")
        with gzip.open(input_path_hg38, "r") as f:
            hg38_dataframe = pd.read_csv(input_path_hg38, sep="\t")


    else:
        print("Nothing found in the file path")

    hg38_dataframe.columns = ["Chromosome", "Start", "End", "mean_coverage"]
    t2t_dataframe.columns = ["Chromosome", "Start", "End", "mean_coverage"]



    return t2t_dataframe, hg38_dataframe

def making_bins(dataframe_t2t, dataframe_hg38):

    mean_coverage_t2t = dataframe_t2t["mean_coverage"].tolist()
    mean_coverage_hg38 = dataframe_hg38["mean_coverage"].tolist()


    coverage_occurences_t2t = Counter(mean_coverage_t2t)
    coverage_occurences_hg38 = Counter(mean_coverage_hg38)
    print(coverage_occurences_hg38)
    print(coverage_occurences_t2t)

    return coverage_occurences_t2t, coverage_occurences_hg38


def make_histogram(occurrences_coverage_t2t, occurrences_coverage_hg38):

    # fig, (ax1, ax2) = plt.subplots(ncols=2, sharey=True)


    # titles = {ax1: "Mean coverage occurrence (CHM13)", ax2: "Mean coverage occurrences (GRCh38)"}
    #
    # ax = plt.gca()

    # ax1.hist(coverage_occurences_t2t, color="blue")
    # ax2.hist(coverage_occurences_hg38, color="red")
    #
    # ax1.set_ylabel("Number of variants in the CHM13")
    # ax2.set_ylabel("Number of variants in the GRCh38")

    plt.hist(occurrences_coverage_t2t, color="blue", bins="auto", label="T2T", alpha=0.5, histtype="step", )
    plt.hist(occurrences_coverage_hg38,color="green", bins="auto", label="GRCh38", alpha=0.25, histtype="stepfilled")
    plt.legend(loc="upper right")
    plt.title("Mean coverage occurences")
    plt.xlabel("Mean coverage")
    plt.ylabel("Number of occurrences")

    plt.savefig("coverage_occurrences_histogram.jpg", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    t2t_dataframe, hg38_dataframe = open_file()
    coverage_occurences_t2t, coverage_occurences_hg38 = making_bins(t2t_dataframe, hg38_dataframe)
    make_histogram(coverage_occurences_t2t, coverage_occurences_hg38)
