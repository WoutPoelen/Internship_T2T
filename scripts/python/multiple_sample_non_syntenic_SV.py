import argparse
import matplotlib
import matplotlib.pyplot as plt

def obtain_structural_variants(args):
    """
    This function obtains the structural variants from the given insertion and deletion files. It then makes the content
    of the files eligible for plotting by removing the whitespaces and splitting the space between the file name and
    number of variants, and lasts adding them to a list of lists.

    :param:
        args.multiple_sample_non_syntenic_INS (txt file): txt file containing per sample the amount of insertions in
        new-syntenic and coding sequences.
        args.multiple_sample_non_syntenic_DEL (txt file): txt file containing per sample the amount of deletions in
        new-syntenic and coding sequences.

    :return:
        list_INS (list): list which contains the amount of insertions per sample.
        list_DEL (list): list which contains the amount of deletions per sample.
    """
    list_INS = []
    list_DEL = []

    # Loops through the file containing the insertions and the file containing the deletions, gets rid of the
    # whitespaces and splits on the space to make a list and then added the list with the sample and the amount of the
    # respective indel to the overall insertion or deletion list
    with open(args.multiple_sample_non_syntenic_INS, "r") as INS_file, \
            open(args.multiple_sample_non_syntenic_DEL, "r") as DEL_file:

        for INS in INS_file:
            split_amount_INS = INS.strip().split(" ")
            list_INS.append(split_amount_INS)

        for DEL in DEL_file:
            split_amount_DEL = DEL.strip().split(" ")
            list_DEL.append(split_amount_DEL)

    return list_INS, list_DEL

def plotting_SV_occurences(list_INS, list_DEL):
    """
    Creates a stacked bar plot showing the amount of non-syntenic insertions and deletions per sample.

    :param:
        list_INS (list): list which contains the amount of insertions per sample.
        list_DEL (list): list which contains the amount of deletions per sample.

    :return:
        Saves the bar plot as Multiple_sample_non_syntenic_SV_occurences.png
    """
    # Makes it possible to be plotted when run on the server
    matplotlib.use("Agg")

    # Sets the size of the figure
    plt.rcParams["figure.figsize"] = [18, 12]

    # Sets the sample names as set keys to compare the deletion and insertion values by looping through the respective
    # indel list and take the value corresponding to the key (sample name)
    keys = [item[0] for item in list_INS]
    del_values = [int(item[1]) for item in list_DEL]
    ins_values = [int(item[1]) for item in list_INS]

    # Sum the values of identical keys to calculate totals for sorting
    total_values = [del_values[i] + ins_values[i] for i in range(len(del_values))]

    # Sort the data by total values by combining the total_values, keys, del values and ins values lists into tuples
    # which contains one element from each list on the same position. Then sort the tuples based on the first element
    # of each tuple (total_values) and then unpack the tuples into separate lists
    sorted_data = sorted(zip(total_values, keys, del_values, ins_values), reverse=True)
    sorted_total, sorted_keys, sorted_del_values, sorted_ins_values = zip(*sorted_data)

    # Creates the stacked bar chart
    plt.bar(sorted_keys, sorted_del_values, label="DEL", color="red")
    plt.bar(sorted_keys, sorted_ins_values, bottom=sorted_del_values, label="INS", color="green")

    # Creates the y-label, change the size and rotation of the sample names, change the size of the numbers and creates
    # the legend and title
    plt.ylabel("Amount of non-syntenic structural variants", fontsize=18)
    plt.xticks(fontsize=12, rotation=80)
    plt.yticks(fontsize=17)
    plt.legend(loc="upper right", fontsize=15)
    plt.title("Amount of non syntenic structural variation per sample", fontsize=20)

    # plt.show(bbox_inches='tight')

    plt.savefig("Multiple_sample_non_syntenic_SV_occurences.png", bbox_inches='tight')

def main(args):
    list_INS, list_DEL = obtain_structural_variants(args)
    plotting_SV_occurences(list_INS, list_DEL)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("multiple_sample_non_syntenic_INS",
                      help="path to the txt file containing the with CDS and non-syntenic regions insertions")
    parser.add_argument("multiple_sample_non_syntenic_DEL",
                        help="path to the txt file containing the with CDS and non-syntenic regions deletions")
    args = parser.parse_args()

    main(args)