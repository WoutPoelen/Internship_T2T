import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

def getting_argument(argument):
    """
    This function recieves two BED files with the amount of times the low coverage regions from their respective
    reference genome overlap with centromeres, transcripts, coding sequences (CDS) and sequential duplications (SD).
    :param:
        argument.intersected_bed_file_T2T_liftover (bed file): bed file containing the T2T to GRCh38 lifted over
        low coverage regions the amount of times it overlaps with a centromere, transcript, coding sequence and
        sequential duplication.

        argument.intersected_bed_file_GRCh38 (bed file): bed file containing the low coverage regions and the amount of
        times it overlaps with a centromere, transcript, coding sequence and sequential duplication.
    :return:
        t2t_liftover_dataframe (dataframe): dataframe containing the low coverage regions who could be lifted over from
        T2T to GRCh38.

        GRCh38_dataframe (dataframe): dataframe containing the low coverage regions from the GRCh38 genome.
    """

    # Open the arguments given bed files and make dataframes from them
    t2t_liftover_dataframe = pd.read_csv(argument.intersected_bed_file_T2T_liftover, sep="\t", encoding="utf-8")
    GRCh38_dataframe = pd.read_csv(argument.intersected_bed_file_GRCh38, sep="\t", encoding="utf-8")

    print("Reading the files")

    # Add column names to the dataframes
    t2t_liftover_dataframe.columns = ["Chromosome","Start", "End", "Coverage", "Category", "Count"]
    GRCh38_dataframe.columns = ["Chromosome", "Start", "End", "Coverage", "Category", "Count"]

    # Turn the amount of times the region overlaps with a difficult category binary. However, some regions would have
    # 7 exomes and that isn't the purpose of this script.
    t2t_liftover_dataframe.loc[t2t_liftover_dataframe["Count"] > 1, "Count"] = 1
    GRCh38_dataframe.loc[GRCh38_dataframe["Count"] > 1, "Count"] = 1

    return t2t_liftover_dataframe, GRCh38_dataframe

def counting_total(liftover_t2t_dataframe, dataframe_GRCh38):
    """
    This function receives the two dataframes from the getting arguments function and gets the regions where nothing
    overlaps. Afterward it gets the total amount of times regions overlap with one of the difficult categories.

    :param:
        liftover_t2t_dataframe (dataframe): dataframe containing the low coverage regions who could be lifted over from
        T2T to GRCh38.

        GRCh38_dataframe (dataframe): dataframe containing the low coverage regions from the GRCh38 genome.

    :return:
        t2t_count_values_dict (dictionary): dictionary containing the total amount of times regions overlap with one (or
        more) of the difficult categories or none of the regions. These regions are lifted over from t2t to grch38.

        grch38_count_values_dict (dictionary): dictionary containing the total amount of times regions overlap with one
        (or more) of the difficult categories or none of the regions. These regions are from the grch38 genome.
    """
    print("Counting overlaps")

    # Groups the dataframe by Start location of the region and the Count to check if the same start location has
    # zero overlap with the categories
    zero_count_per_start_t2t = liftover_t2t_dataframe.groupby('Start')['Count'].apply(lambda x: (x == 0).all())
    zero_count_per_start_grch38 = dataframe_GRCh38.groupby('Start')['Count'].apply(lambda x: (x == 0).all())

    # Get the 'Start' positions where all counts are zero and turns it into a list
    start_positions_with_zero_counts_t2t = zero_count_per_start_t2t[zero_count_per_start_t2t].index.tolist()
    start_positions_with_zero_counts_grch38 = zero_count_per_start_grch38[zero_count_per_start_grch38].index.tolist()

    # The amount of not overlapping regions
    total_no_overlaps_T2T = len(start_positions_with_zero_counts_t2t)
    total_no_overlaps_GRCh38 = len(start_positions_with_zero_counts_grch38)

    # Group the DataFrame by "Category" and sum the "Count" column to get the amount of regions that overlap with the
    # categories
    t2t_count_values = liftover_t2t_dataframe.groupby("Category")["Count"].sum()
    grch38_count_values = dataframe_GRCh38.groupby("Category")["Count"].sum()

    # Converts the dataframes into dictionaries a
    t2t_count_values_dict = t2t_count_values.to_dict()
    grch38_count_values_dict = grch38_count_values.to_dict()

    # Adds the amount of non overlapping regions the dictionaries
    t2t_count_values_dict["no overlap"] = total_no_overlaps_T2T
    grch38_count_values_dict["no overlap"] = total_no_overlaps_GRCh38

    # Hardcoded failed liftover value. Change this to the appropriate value
    # t2t_count_values_dict["Failed liftover"] = 207089
    # grch38_count_values_dict["Failed liftover"] = 0

    return t2t_count_values_dict, grch38_count_values_dict

def making_barplot(t2t_count_values, grch38_count_values):
    """
    This functions recieves the dictionaries from counting total function and combines them into one dictionary.
    Then makes a list with the category names and gets the values in list form from their respective dictionary and
    reference genome. The lists with keys and values get plotted into a barplot.
    :param:
        t2t_count_values_dict (dictionary): dictionary containing the total amount of times regions overlap with one (or
        more) of the difficult categories or none of the regions. These regions are lifted over from t2t to grch38.
        grch38_count_values_dict (dictionary): dictionary containing the total amount of times regions overlap with one
        (or more) of the difficult categories or none of the regions. These regions are from the grch38 genome.
    :return:
        The plot is saved in low_coverage_categories_barplot.png
    """
    print("Generating barplot")
    
    # Combines the two dictionaries to make it easier to plot the grouped barplot
    complete_dictionary = {"T2T": t2t_count_values, "GRCh38": grch38_count_values}

    # Makes a list of the categories to use as ticks for the plot
    categories = list(complete_dictionary["T2T"].keys())

    # Makes a list of the values of the dictionaries by going to the T2T key and looping through the categories and
    # getting their values
    t2t_values_list = [complete_dictionary["T2T"][category] for category in categories]
    grch38_values_list = [complete_dictionary["GRCh38"][category] for category in categories]

    # Give an x-value to the bar plots, so the bar plots won't stack on top of each other when added or subtracted from
    x = np.arange(len(categories))

    # Makes matplotlib work on the server.
    matplotlib.use("Agg")

    # Makes subplots
    fig, ax = plt.subplots(figsize=(8, 6))

    # Generates the barplots with the x-value depending on the reference genome
    # SD is first so its x-value is 0 and for t2t it's -0.2 and grch38 is 0.2.
    ax.bar(x - 0.2, t2t_values_list, width=0.4, label = "T2T", color = "blue")
    ax.bar(x + 0.2, grch38_values_list, width=0.4, label = "GRCh38", color = "green")

    # Set the locations of the category names
    ax.set_xticks(x)

    # If a limit on the y_axis is necessary. Remove the comment
    # ax.set_ylim(0, 180000)

    # Set the categories as labels beneath the bars
    ax.set_xticklabels(categories)

    # Generate the legend, title and x/y label
    ax.legend()
    ax.set_title("Amount of low coverage regions overlapping with difficult regions")
    ax.set_ylabel("Amount of overlapping regions")
    ax.set_xlabel("Categories")

    # Saves the figure as the following png file
    plt.savefig("low_coverage_categories_barplot.png")

    # Shows the plot. Get rid of comment if run in pycharm
    # plt.show()


def main(args):
    t2t_liftover_dataframe, GRCh38_dataframe = getting_argument(args)
    t2t_count_values_dict, grch38_count_values_dict = (counting_total(t2t_liftover_dataframe, GRCh38_dataframe))
    making_barplot(t2t_count_values_dict, grch38_count_values_dict)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Gets the bed files containing the low coverage regions and the "
                                                 "difficult categories (structural duplications, transcription, "
                                                 "coding sequences and centromeres) and counts the total amount of "
                                                 "times regions overlap with the categories. "
                                                 "Plots those amounts in a bar plot")
    parser.add_argument("intersected_bed_file_T2T_liftover",
                        help="Path to the liftover (t2t to GRCh38) bed file containing the amount of categories a "
                             "region overlaps with",
                        metavar="Intersected bed file from the T2T to GRCh38 liftover")
    parser.add_argument("intersected_bed_file_GRCh38",
                        help="Path to the GRCh38 bed file containing the amount of categories a region overlaps with",
                        metavar="Intersected bed file from the GRCh38 file")
    args = parser.parse_args()
    main(args)