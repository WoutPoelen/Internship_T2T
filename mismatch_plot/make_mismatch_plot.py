import os
import bamnostic as bs  #bamnostic-1.1.10
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# C:\Users\Z659216\Desktop\BAM_files\P3-D10.aligned.HiFi.bam
# C:\Users\Z659216\Desktop\BAM_files\P3-D10.haplotagged.bam

def read_files():
    """This function asks for two paths: one to the t2t bam file and one to the hg38 file of the same sample and
    then adds the reads to list depending on which file the sample is from and then returns the lists.

    Args:
        No arguments

    Returns:
        hg38_reads(list): list of reads from the hg38 bam file.
        t2t_reads(list): list of reads from the t2t bam file.
    """

    t2t_reads = []
    hg38_reads = []
    input_path_t2t_bam = input("Give the path to the t2t BAM file: ")
    print("given path: " + input_path_t2t_bam)

    if os.path.exists(input_path_t2t_bam):
        print("File exists")
        print("Processing the reads")
        print(" ")
        # Opens the file and turns it into a pandas dataframe.
        with bs.AlignmentFile(input_path_t2t_bam, "rb") as bam_t2t:
            for read in bam_t2t.head(n=15000):
                t2t_reads.append(read)

        bam_t2t.close()


    # Prints the string if the file path doesn't exist
    else:
        print("The specified path does NOT exist!")

    input_path_hg38_bam = input("Give the path to the GRCh38 BAM file: ")
    print("given path: " + input_path_hg38_bam)

    if os.path.exists(input_path_hg38_bam):
        print("File exists")
        print("Processing the reads")
        print(" ")
        with bs.AlignmentFile(input_path_hg38_bam, "rb") as bam_hg38:
            for read in bam_hg38.head(n=15000):
                hg38_reads.append(read)
        bam_hg38.close()

    else:
        print("The specified path does NOT exist!")



    return hg38_reads, t2t_reads

def get_necessary_data(reads_hg38, reads_t2t):
    """This function gets the lists with reads and searches per read for the NM tag (contains the number of mismatches)
    and the length of the alignment part of the read. Then it calculates the mismatch rate by dividing the number of
    mismatches by the length of the alignment part of the read. Returns the lists with the mismatches rates for the
    two files.

    Args:
        reads_hg38 (list): list of reads from the hg38 bam file.
        reads_t2t (list): list of reads from the t2t bam file.

    Returns:
        mismatch_rate_list_hg38(list): list of mismatches rates from the hg38 bam file.
        mismatch_rate_list_t2t(list): list of mismatches rates from the t2t bam file.
    """

    mismatch_rate_list_t2t = []
    mismatch_rate_list_hg38 = []

    for read in reads_hg38:
        mismatch_tag = read.get_tag("NM")
        if mismatch_tag is not None:
            mismatch_rate = mismatch_tag / read.query_alignment_length
            mismatch_rate_list_hg38.append(mismatch_rate)

    for read in reads_t2t:
        mismatch_tag = read.get_tag("NM")
        if mismatch_tag is not None:
            mismatch_rate = mismatch_tag / read.query_alignment_length
            mismatch_rate_list_t2t.append(mismatch_rate)
    mean_mismatch_hg38 = sum(mismatch_rate_list_hg38) / len(reads_hg38)
    mean_mismatch_t2t = sum(mismatch_rate_list_t2t) / len(reads_t2t)


    return mean_mismatch_hg38, mean_mismatch_t2t

def plot_boxplot_comparison(mismatch_mean_hg38, mismatch_mean_t2t):
    """Uses the lists with mismatch rate to plot a violin plot. Use is for comparison between the two reference genomes.

    Args:
        mismatch_rate_t2t_list(list): list of mismatches rates from the t2t bam file.
        mismatch_rate_hg38_list(list): list of mismatches rates from the hg38 bam file.
    """
    dataframe_combined_mismatches = pd.DataFrame({"T2T": [mismatch_mean_t2t], "Hg38": [mismatch_mean_hg38]},
                                                 columns=["T2T", "Hg38"])
    fig, ax =plt.subplots(figsize=(5, 4))

    ax.scatter(x=dataframe_combined_mismatches.keys(), y=dataframe_combined_mismatches.values)
    ax.set_xlim(-0.50, 1.50)
    ax.set_ylim(0, 0.04)
    ax.set_xticklabels(["T2T", "Hg38"])

    plt.title("Average mismatch rate")
    plt.savefig("mismatch_rate_comparison.jpg")

    plt.show()


if __name__ == "__main__":
    hg38_reads, t2t_reads = read_files()
    mean_mismatch_hg38, mean_mismatch_t2t = get_necessary_data(hg38_reads, t2t_reads)
    plot_boxplot_comparison(mean_mismatch_hg38, mean_mismatch_t2t)