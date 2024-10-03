import argparse
import os
import bamnostic as bs
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

def read_files(arguments):
    """
    This function asks for two paths: one to the t2t bam file and one to the hg38 file of the same sample and
    then adds the reads to list depending on which file the sample is from and then returns the lists.

    Args:
        No arguments

    Returns:
        hg38_reads(list): list of reads from the hg38 bam file.
        t2t_reads(list): list of reads from the t2t bam file.
    """
    # lists where the reads will be appended into
    t2t_reads = []
    hg38_reads = []

    print("Reading the files")
    #     # Opens the file and loops through it so it can add the reads into the bam_t2t list
    with bs.AlignmentFile(arguments.T2T_BAM_file, "rb") as bam_t2t:
        for read in bam_t2t.head(10):
            t2t_reads.append(read)

        bam_t2t.close()

    #     # Opens the file and loops through it so it can add the reads into the bam_hg38 list
        with bs.AlignmentFile(arguments.GRCh38_BAM_file, "rb") as bam_hg38:
            for read in bam_hg38.head(10):
                hg38_reads.append(read)

        bam_hg38.close()


    return hg38_reads, t2t_reads

def get_necessary_data(reads_hg38, reads_t2t):
    """
    This function gets the lists with reads and searches per read for the NM tag (contains the number of mismatches)
    and the length of the alignment part of the read. Then it calculates the mismatch rate by dividing the number of
    mismatches by the length of the alignment part of the read. Returns the lists with the mismatches rates for the
    two files.

    Args:
        reads_hg38 (list): list of reads from the hg38 bam file.
        reads_t2t (list): list of reads from the t2t bam file.

    Returns:
        mean_mismatch_hg38(float): the average of the mismatch rate across the entire file. Calculated by dividing the
        sum of mismatches by the length of the list with reads.
        mean_mismatch_t2t(float): the average of the mismatch rate across the entire file. Calculated by dividing the
        sum of mismatches by the length of the list with reads
    """

    # The lists that will contain the mismatch rates
    mismatch_rate_list_t2t = []
    mismatch_rate_list_hg38 = []

    print("Calculating the mismatch rates")

    # Loops through the list with the reads from the hg38 BAM file and gets the number of mismatches from the reads
    # by getting the NM tag. Calculates the mismatch rate of the read by dividing the number of mismatches in the read
    # by the length of the aligned portion of the read. Only if read has a tag
    for read in reads_hg38:
        mismatch_tag = read.get_tag("NM")
        if mismatch_tag is not None:
            mismatch_rate = mismatch_tag / read.query_alignment_length
            mismatch_rate_list_hg38.append(mismatch_rate)

    # Loops through the list with the reads from the t2t BAM file and gets the number of mismatches from the reads
    # by getting the NM tag. Calculates the mismatch rate of the read by dividing the number of mismatches in the read
    # by the length of the aligned portion of the read. Only if read has a tag
    for read in reads_t2t:
        mismatch_tag = read.get_tag("NM")
        if mismatch_tag is not None:
            mismatch_rate = mismatch_tag / read.query_alignment_length
            mismatch_rate_list_t2t.append(mismatch_rate)

    # Calculates the average mismatch rate of the files by dividing the sum of mismatches in their respective list by
    # the length of the list with the reads from that file.
    mean_mismatch_hg38 = sum(mismatch_rate_list_hg38) / len(reads_hg38)
    mean_mismatch_t2t = sum(mismatch_rate_list_t2t) / len(reads_t2t)

    return mean_mismatch_hg38, mean_mismatch_t2t

def plot_boxplot_comparison(mismatch_mean_hg38, mismatch_mean_t2t):
    """
    Uses the lists with mismatch rate to plot a violin plot. Use is for comparison between the two reference genomes.

    Args:
        mismatch_mean_hg38(float): the average of the mismatch rate across the entire file. Calculated by dividing the
        sum of mismatches by the length of the list with reads.
        mismatch_mean_t2t(float): the average of the mismatch rate across the entire file. Calculated by dividing the
        sum of mismatches by the length of the list with reads.

    Returns:
        saves the plot as mismatch_rate_comparison.jpg.
    """

    matplotlib.use("Agg")

    print("Generating the plot")

    colors = ["blue", "red"]
    # Makes a dataframe containing the average mismatches of the files.
    dataframe_combined_mismatches = pd.DataFrame = pd.DataFrame({"T2T": [mismatch_mean_t2t],
                                                                 "Hg38": [mismatch_mean_hg38]
                                                                 })


    # Sets the figure up with the width and height
    fig, ax =plt.subplots(figsize=(5, 4))

    # Creates a scatter plot with the name of the reference genomes in the x-axis and the average mismatch rates in
    # the y-axis in the figure
    ax.bar(x=dataframe_combined_mismatches.columns, height=dataframe_combined_mismatches.iloc[0], color=colors,
           width=0.5)

    # Sets a limit to the x-axis and y-axis
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(0, 0.04)

    # Sets the x labels up
    ax.set_xticklabels(["T2T", "Hg38"])

    # Sets the title up
    plt.title("Average mismatch rate")

    # Saves the figure in a jpg file
    plt.savefig("mismatch_rate_comparison.png")

    print("The barplot has been generated and is saved as mismatch_rate_comparison.png")

    # Makes sure the plot is shown. Commented out because file wouldn't be saved on the server. Get rid of the comment
    # when you want to see the file outside the server
    # plt.show()

def main(args):
    hg38_reads, t2t_reads = read_files(args)
    mean_mismatch_hg38, mean_mismatch_t2t = get_necessary_data(hg38_reads, t2t_reads)
    plot_boxplot_comparison(mean_mismatch_hg38, mean_mismatch_t2t)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("T2T_BAM_file",
                        help="Path to the BAM file of the sample that is mapped to T2T",
                        metavar="the T2T input BAM file")
    parser.add_argument("GRCh38_BAM_file",
                        help="Path to the BAM file of the sample that is mapped to T2T",
                        metavar="the GRCh38 input BAM file")
    args = parser.parse_args()

    main(args)