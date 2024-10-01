import pandas as pd
import os

from coverage_occurrence_histogram.make_coverage_mean_histogram import open_file


def openfile():
    # Asks for the input path for the t2t bed regions file
    input_path_t2t = input("Give the path to the t2t bed regions file: ")
    print("given path: " + input_path_t2t)

    # Checks if there is a file in the file path
    if os.path.exists(input_path_t2t):
        print("File exists")
        print(" ")
        # Opens the file with gunzip and turns it into a pandas dataframe
        t2t_dataframe = pd.read_csv(input_path_t2t, sep="\t", compression="gzip")

        # Asks for the input path for the t2t BAM file
    input_path_hg38 = input("Give the path to the GRCh38 bed regions file: ")
    print("given path: " + input_path_hg38)

    # Checks if there is a file in the file path
    if os.path.exists(input_path_hg38):
        print("File exists")
        print(" ")
        # Opens the file with gunzip and turns it into a dataframe
        hg38_dataframe = pd.read_csv(input_path_hg38, sep="\t", compression="gzip")

    return t2t_dataframe, hg38_dataframe

def process_dataframes(dataframe_t2t, dataframe_hg38):
    dataframe_t2t.columns = ["Chromosome", "Start", "End", "mean_coverage"]
    dataframe_hg38.columns = ["Chromosome", "Start", "End", "mean_coverage"]

    hg38_dataframe = dataframe_hg38[(dataframe_hg38["End"] <= 5000000) & (dataframe_hg38["Chromosome"] == "chr1") &
                                    (dataframe_hg38["mean_coverage"] <= 10)]
    t2t_dataframe = dataframe_t2t[(dataframe_t2t["End"] <= 5000000) & (dataframe_t2t["Chromosome"] == "chr1") &
                                  (dataframe_t2t["mean_coverage"] <= 10)]

    return t2t_dataframe, hg38_dataframe

def write_to_file(t2t_dataframe, hg38_dataframe):

    t2t_dataframe.to_csv("t2t_low_coverage.bed", sep="\t", index=False, header=False)
    hg38_dataframe.to_csv("hg38_low_coverage.bed", sep="\t", index=False, header=False)


def main():
    t2t_dataframe, hg38_dataframe = open_file()
    t2t_dataframe, hg38_dataframe = process_dataframes(t2t_dataframe, hg38_dataframe)
    write_to_file(t2t_dataframe, hg38_dataframe)
if __name__ == "__main__":
    main()