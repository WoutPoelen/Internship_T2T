import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import gzip
import os

def openfile():
    # Asks for the input path for the t2t bed regions file
    input_path_t2t = input("Give the path to the t2t bed regions file: ")
    print("given path: " + input_path_t2t)

    # Checks if there is a file in the file path
    if os.path.exists(input_path_t2t):
        print("File exists")
        print(" ")
        # Opens the file with gunzip and turns it into a pandas dataframe
        t2t_dataframe = pd.read_csv(input_path_t2t, sep="\t")

        # Asks for the input path for the t2t BAM file
    input_path_hg38 = input("Give the path to the GRCh38 bed regions file: ")
    print("given path: " + input_path_hg38)

    # Checks if there is a file in the file path
    if os.path.exists(input_path_hg38):
        print("File exists")
        print(" ")
        # Opens the file with gunzip and turns it into a dataframe
        hg38_dataframe = pd.read_csv(input_path_hg38, sep="\t")

    t2t_dataframe.columns = ["Chromosome", "Start", "End", "mean_coverage"]
    hg38_dataframe.columns = ["Chromosome", "Start", "End", "mean_coverage"]

    hg38_dataframe = hg38_dataframe[(hg38_dataframe["End"] <= 5000000) & (hg38_dataframe["Chromosome"] == "chr1") &
                                    (hg38_dataframe["mean_coverage"] <= 10)]
    t2t_dataframe = t2t_dataframe[(t2t_dataframe["End"] <= 5000000) & (t2t_dataframe["Chromosome"] == "chr1") &
                                  (t2t_dataframe["mean_coverage"] <= 10)]

    print(t2t_dataframe)
    print(hg38_dataframe)
    return t2t_dataframe, hg38_dataframe

def plot_low_coverage(dataframe_t2t, dataframe_hg38):
   dataframe_t2t["Reference genome"] = "T2T"
   dataframe_hg38["Reference genome"] = "GRCh38"


   figs, ax = plt.subplots(figsize=(12, 5))
   plt.scatter(dataframe_t2t["Start"], dataframe_t2t["Reference genome"], s=1)
   plt.scatter(dataframe_hg38["Start"], dataframe_hg38["Reference genome"], s=1)
   # plt.vlines(x=dataframe_t2t["Start"], color="r", ymin=0, ymax=1)
   # plt.vlines(x=dataframe_hg38["Start"], color="b", ymin=0, ymax=1)

   # x_ticks = np.arange(0, dataframe_t2t["Start"].max(), step=10000)  # Adjust the step size
   # ax.set_xticks(x_ticks)

   plt.show()


def main():
    t2t_dataframe, hg38_dataframe = openfile()
    plot_low_coverage(t2t_dataframe, hg38_dataframe)


if __name__ == "__main__":
    main()