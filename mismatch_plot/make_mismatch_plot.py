import os
import bamnostic as bs  #bamnostic-1.1.10
import pandas as pd
import matplotlib.pyplot as plt

# C:\Users\Z659216\Desktop\BAM_files\P3-D10.aligned.HiFi.bam
# C:\Users\Z659216\Desktop\BAM_files\P3-D10.haplotagged.bam

def read_files():
    """read the input files"""
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
            for read in bam_t2t.head(n=20000):
                t2t_reads.append(read)
        bam_t2t.close()


    # Prints the string if the file path doesn't exist
    else:
        print("The specified path does NOT exist!")

    input_path_hg38_bam = input("Give the path to the hg38 BAM file: ")
    print("given path: " + input_path_hg38_bam)

    if os.path.exists(input_path_hg38_bam):
        print("File exists")
        print("Processing the reads")
        print(" ")
        with bs.AlignmentFile(input_path_hg38_bam, "rb") as bam_hg38:
            for read in bam_hg38.head(n=20000):
                hg38_reads.append(read)
        bam_hg38.close()

    else:
        print("The specified path does NOT exist!")



    return hg38_reads, t2t_reads

def get_necessary_data(reads_hg38, reads_t2t):
    """gets the alignment file and extracts the alignment length and the number of mismatches"""
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

    return mismatch_rate_list_hg38, mismatch_rate_list_t2t

def plot_boxplot_comparison(mismatch_rate_hg38_list, mismatch_rate_t2t_list):
    """gets the alignment length and the number of mismatches"""

    df_mismatches = pd.DataFrame({"T2T": mismatch_rate_t2t_list, "HG38": mismatch_rate_hg38_list})
    fig, ax = plt.subplots()
    ax.boxplot(df_mismatches, tick_labels=df_mismatches.keys())
    plt.show()

if __name__ == "__main__":
    hg38_reads, t2t_reads = read_files()
    mismatch_rate_list_hg38, mismatch_rate_list_t2t =get_necessary_data(hg38_reads, t2t_reads)
    plot_boxplot_comparison(mismatch_rate_list_hg38, mismatch_rate_list_t2t)