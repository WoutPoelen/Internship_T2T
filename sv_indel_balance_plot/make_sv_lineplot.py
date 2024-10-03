import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

def open_sv_files(arguments):
    """This function reads the t2t sv file and turns it into a dataframe. Checks if the row is an insertion of deletion
     depending on the length of the variation and makes the length of the deletions absolute.

     :param:
        argument.T2T_filtered_txt_file (txt file): text file with the Structural Variants (sv) from the sample mapped
        to the T2T reference genome.
        argument.GRCh38_filtered_txt_file (txt file): text file with the Structural Variants (sv) from the sample mapped
        to the GRCh38 (Hg38) reference genome.

     :return:
        df_t2t_insertions (pandas dataframe): dataframe containing insertions from the t2t file.
        df_t2t_deletions_absolute (pandas dataframe): dataframe containing deletions from the t2t file in absolute
        format.
        df_hg38_insertions (pandas dataframe): dataframe containing insertions from the hg38 file.
        df_hg38_deletions_absolute (pandas dataframe): dataframe containing deletions from the hg38 file in absolute
        format.
    """
    # Opens the argument files and makes them into dataframes
    t2t_dataframe = pd.read_csv(arguments.T2T_filtered_txt_file, sep=" ", encoding="utf-8")
    hg38_dataframe = pd.read_csv(arguments.GRCh38_filtered_txt_file, sep=" ", encoding="utf-8")

    # Gives the column names to the dataframe
    t2t_dataframe.columns = ["CHROM", "POS", "FILTER", "GT", "LENGTH"]
    hg38_dataframe.columns = ["CHROM", "POS", "FILTER", "GT", "LENGTH"]

    # Separates the insertions and deletions by looking at the lengths if the indels (positive length means insertion
    # negative length means deletions)
    df_t2t_insertions = t2t_dataframe[t2t_dataframe["LENGTH"] > 0]
    df_t2t_deletions = t2t_dataframe[t2t_dataframe["LENGTH"] < 0]

    df_hg38_insertions = hg38_dataframe[hg38_dataframe["LENGTH"] > 0]
    df_hg38_deletions = hg38_dataframe[hg38_dataframe["LENGTH"] < 0]

    # Turns the length of the deletions absolute
    df_t2t_deletions_absolute = df_t2t_deletions["LENGTH"].abs()
    df_hg38_deletions_absolute = df_hg38_deletions["LENGTH"].abs()

    return df_t2t_insertions, df_t2t_deletions_absolute, df_hg38_insertions, df_hg38_deletions_absolute


def get_and_compare_lengths(df_t2t_insertions, t2t_df_deletions_absolute, df_hg38_insertions,
                            df_hg38_deletions_absolute):
    """
    This functions gets the dataframes with insertions or deletions from both the t2t and hg38 reference genome.
    It then goes through the dataframes and gets the amount of variations whose length is in between certain lengths.

    :param:
        df_t2t_insertions (dataframe): dataframe containing insertions from the t2t file.
        t2t_df_deletions_absolute (dataframe): dataframe containing deletions whose length has been made absolute
        from the t2t file.
        df_hg38_insertions (dataframe): dataframe containing insertions from the hg38 file.
        hg38_df_deletions (dataframe): dataframe containing deletions from the hg38 file.

    Returns:
        t2t_ins_distance (list): list with the amount of insertions from the t2t file whose length is in between certain
        lengths.
        t2t_del_distances (list): list with the amount of deletions from the t2t file whose length is in between certain
        lengths.
        hg38_ins_distances (list): list with the amount of insertions from the hg38 file whose length is in between
        certain lengths.
        hg38_del_distances (list): list with the amount of deletions from the hg38 file whose length is in between
        certain lengths.
        highest_count (int): the highest count of variations in compared to the groups.
    """

    # Creates empty lists that will be used for the lengths of the into group separated variations.
    hg38_ins_distances = []
    hg38_del_distances = []
    t2t_ins_distances = []
    t2t_del_distances = []

    # Separates the insertions from the hg38 file based on the length of the insertions and adds the length of the
    # subset into the list that will be used as values for the plot. So hg38_ins_distances contains at the end 12
    # numbers corresponding with the amount of variations that have a certain length.
    hg38_ins_distances.extend(
        [len(df_hg38_insertions["LENGTH"][df_hg38_insertions["LENGTH"].between(start, end)]) for start, end in
         [(30, 50), (50, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 750), (750, 1000), (1000, 2000),
          (2000, 5000), (5000, 10000)] + [(10000, float('inf'))]])

    # Separates the insertions from the hg38 file based on the length of the deletions and adds the length of the
    # subset into the list that will be used as values for the plot. So hg38_del_distances contains at the end 12
    # numbers corresponding with the amount of variations that have a certain length.
    hg38_del_distances.extend(
        [len(df_hg38_deletions_absolute[df_hg38_deletions_absolute.between(start, end)]) for start,
        end in [(30, 50), (50, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 750), (750, 1000),
                (1000, 2000), (2000, 5000), (5000, 10000)] + [(10000, float('inf'))]])
    # ----------------------------------------------------------------------------------------------------

    # Separates the insertions from the t2t file based on the length of the insertions and adds the length of the
    # subset into the list that will be used as values for the plot. So hg38_ins_distances contains at the end 12
    # numbers corresponding with the amount of variations that have a certain length.

    t2t_ins_distances.extend(
        [len(df_t2t_insertions["LENGTH"][df_t2t_insertions["LENGTH"].between(start, end)]) for start, end in
         [(30, 50), (50, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 750), (750, 1000), (1000, 2000),
          (2000, 5000), (5000, 10000)] + [(10000, float('inf'))]])

    # Separates the insertions from the t2t file based on the length of the deletions and adds the length of the
    # subset into the list that will be used as values for the plot. So hg38_ins_distances contains at the end 12
    # numbers corresponding with the amount of variations that have a certain length.
    t2t_del_distances.extend(
        [len(t2t_df_deletions_absolute[t2t_df_deletions_absolute.between(start, end)]) for start, end in
         [(30, 50), (50, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 750), (750, 1000), (1000, 2000),
          (2000, 5000), (5000, 10000)] + [(10000, float('inf'))]])

    # Gets the highest number of variants in a group from both the insertion and deletion dataframe. Then gets the
    # highest number between the two. T2T is the only one to be checked, because it has more reads on average.
    # This is used for the y-axis of the plots.
    max_deletions = max(t2t_del_distances)
    max_insertions = max(t2t_ins_distances)
    highest_count = max(max_deletions, max_insertions)

    return t2t_ins_distances, t2t_del_distances, hg38_ins_distances, hg38_del_distances, highest_count


def plot(t2t_ins_distances, t2t_del_distances, hg38_ins_distances, hg38_del_distances, highest_count):
    """This function plots one figure containing two different plots. The first plot has the indels balance with the
    CHM13 reference genome and the second plot has the indel balance with the Hg38 reference genome.

    Args:
        t2t_ins_distances (list): list of the lengths of the insertions found with the t2t reference genome.
        t2t_del_distances (list): list of the lengths of the deletions found with the t2t reference genome.
        hg38_ins_distances (list): list of the lengths of the insertions found with the hg38 reference genome.
        hg38_del_distances (list): list of the lengths of the deletions found with the hg38 reference genome.

    Returns:
        nothing
    """

    # Makes the plotting work on the server
    matplotlib.use("Agg")

    # List with the different groups of indel lengths. Used for the x-axis.
    lengths = ["30-50", "50-100 ", "100-200", "200-300", "300-400", "400-500", "500-750", "750-1k", "1k-2k", "2k-5k",
               "5k-10k", "10k+"]

    # Sets the font size for the plots.
    # plt.rcParams.update({'font.size': 5.75})

    # Indicates that two plots (ax1, ax2) will be made, that they will be next to each other (2 columns) and that
    # they share the same y-axis.
    fig, (ax1, ax2) = plt.subplots(ncols=2, sharey=True, figsize=(10,8))

    # The titles for the two plots. In a dictionary that is later given to a for loop that sets everything up.
    titles = {ax1: "SV indel balance (CHM13)", ax2: "SV indel balance (GRCh38)"}

    # To return the current active axes. Used to set extra things for the plots up.
    ax = plt.gca()

    # Setting up the first plot. This has the lengths list as the x-axis, amount of insertions and deletions from the
    # t2t file separate on the y-axis, marker to get a dot when the x and y-axis cross, color to dictate the color and
    # label to get that string in the legend for that plot (line in this case.)
    ax1.plot(lengths, t2t_ins_distances, marker= "o", color="blue", label="INS")
    ax1.plot(lengths, t2t_del_distances, marker="o", color="red", label="DEL")

    # Setting up the second plot. This has the lengths list as the x-axis, amount of insertions and deletions from the
    # Hg38 file separate on the y-axis, marker to get a dot when the x and y-axis cross, color to dictate the color and
    # label to get that string in the legend for that plot (line in this case).
    ax2.plot(lengths, hg38_ins_distances, marker= "o", color="blue", label="INS")
    ax2.plot(lengths, hg38_del_distances, marker= "o", color="red",label="DEL")

    # Set the y-labels for both of the plots.
    ax1.set_ylabel("Number of variants in the CHM13")
    ax2.set_ylabel("Number of variants in the GRCh38")

    # Setting the limit of the y_axis. Use the highest count from the previous function to have the limit change with
    # each dataset.
    ax.set_ylim([0, highest_count + 1000])

    # Loops through both of the plots and adds the following things the plots: the title corresponding with the plot
    # from the dictionary, a legend, a grid, rotated the x-values 40 degrees and the x-label
    for ax in [ax1, ax2]:
        ax.set_title(titles[ax])
        ax.legend(loc="upper right")
        ax.grid(True)
        ax.tick_params(labelrotation=40, axis="x")
        ax.set_xlabel("Length")

    # Saves the plot as the given file.
    plt.savefig("SV_indel_comparison.png", bbox_inches="tight")

    print("The plot is successfully generated and saved in SV_indel_comparison.png")
    # Makes sure the plot is shown. Commented out because file wouldn't be saved on the server. Get rid of the comment
    # when you want to see the file outside the server
    # plt.show()

def main(args):
    df_t2t_insertions, t2t_df_deletions_absolute, df_hg38_insertions, df_hg38_deletions_absolute = open_sv_files(args)
    t2t_ins_distances, t2t_del_distances, hg38_ins_distances, hg38_del_distances, highest_count = (
        get_and_compare_lengths(df_t2t_insertions, t2t_df_deletions_absolute, df_hg38_insertions,
                                df_hg38_deletions_absolute))
    plot(t2t_ins_distances, t2t_del_distances, hg38_ins_distances, hg38_del_distances, highest_count)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("T2T_filtered_txt_file",
                        help="Path to the filtered BED file from process_regions_file.py containing "
                             "the regions of 500 bp and the average coverage lower than 10 from the T2T file")
    parser.add_argument("GRCh38_filtered_txt_file",
                        help="Path to the filtered BED file from process_regions_file.py containing the regions "
                             "of 500 bp and the average coverage lower than 10 from the GRCh38 file")
    args = parser.parse_args()
    main(args)