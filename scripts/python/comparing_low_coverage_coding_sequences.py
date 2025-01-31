import argparse
import matplotlib
import numpy
import pandas as pd
from matplotlib import pyplot as plt

def get_files(arguments):
    """
    This function receives two Browser Extensible Data files (BED) containing the genes whose Coding Sequences (CDS)
    overlap with low coverage regions from their respective reference genome and puts the Chromosome, Start, End and
    gene_id of each overlap in a dataframe. It also adds _Y to gene_ids on the Y chromosomes and drops
    duplicate Gene_ids.

    :param:
        arguments.T2T_intersected_CDS_Only (BED file): BED file containing the genes in T2T RefSeq
        whose Coding Sequences (CDS) overlap with low coverage region in the T2T Categorical CDS only file.
        arguments.GRCh38_intersected_CDS_only (BED file): BED file containing the genes in GRCh38 (hg38) RefSeq
        whose CDS overlap with low coverage region in the GRCh38 Categorical CDS Only file.

    :return:
        t2t_dataframe (dataframe): dataframe containing the chromosome, start, end and Gene_id of the CDS regions in the
        T2T intersect file.
        hg38_dataframe: dataframe containing the chromosome, start, end and Gene_id of the CDS regions in the
        GRCh38 intersect file.
    """
    print("Reading files", "\n")

    column_names = ["Chromosome", "Start", "End", "Gene_id"]

    # Turns the two bed files into pandas dataframes
    t2t_dataframe = pd.read_csv(arguments.T2T_intersected_CDS_Only, sep="\t", encoding="utf-8", names=column_names)
    hg38_dataframe = pd.read_csv(arguments.GRCh38_intersected_CDS_only, sep="\t", encoding="utf-8", names=column_names)

    # Add _Y to the end of the genes on the Y chromosome, because there are genes that are on both X and Y chromosome,
    # but only on the Y chromosome is it low coverage
    t2t_condition = (t2t_dataframe["Chromosome"] == "chrY")
    hg38_condition = (hg38_dataframe["Chromosome"] == "chrY")
    t2t_dataframe.loc[t2t_condition, "Gene_id"] = t2t_dataframe.loc[t2t_condition, "Gene_id"].astype(str) + "_Y"
    hg38_dataframe.loc[hg38_condition, "Gene_id"] = hg38_dataframe.loc[hg38_condition, "Gene_id"].astype(str) + "_Y"

    # Splits the combined geneIDs that were made in step four in low_coverage_coding_sequence_regions.sh through the bedtools merge.
    # Making the geneIDs a list containing the geneIDS in that region
    t2t_dataframe["Gene_id"] = t2t_dataframe["Gene_id"].str.split(",")
    hg38_dataframe["Gene_id"] = hg38_dataframe["Gene_id"].str.split(",")

    # Makes separate rows for the genes in the lists with multiple geneIDS. This may cause the amount of unique genes to
    # decrease (if there are multiple combination with the same gene: IGLJ1,IGLL5 and IGLC1,IGLL5)
    t2t_dataframe = t2t_dataframe.explode("Gene_id")
    hg38_dataframe = hg38_dataframe.explode("Gene_id")

    return t2t_dataframe, hg38_dataframe


def sort_files(t2t_dataframe, hg38_dataframe):
    """
    This function takes the dataframes made in the previous function and compares the gene_ids between the dataframes to
    see which gene_ids have a CDS which overlaps with a low coverage region.

    Adds the gene_id based on the comparison result:
        - If the gene_id is in both dataframes, it is added to common_low_coverage_regions_genes_length.
        - If the gene_id is only in the T2T dataframe, it is added to t2t_unique_low_coverage_regions_length.
        - If the gene_id is only in the GRCh38 dataframe, it is added to GRCh38_unique_low_coverage_regions_genes_length.

    The lengths of common_low_coverage_regions_genes_length, t2t_unique_low_coverage_regions_length, and
    GRCh38_unique_low_coverage_regions_genes_length are then calculated and saved.

    :param:
        t2t_dataframe (dataframe): Dataframe containing the chromosome, start, end and Gene_id of the CDS regions in
        T2T with every Gene_id only appearing once.
        hg38_dataframe: Dataframe containing the chromosome, start, end and Gene_id of the CDS regions in
        GRCh38 with every Gene_id only appearing once.

    :return:
        common_low_coverage_regions_genes_length (int): the amount of genes who have a low coverage coding sequence in
        both reference genomes. Used for the barplot
        t2t_unique_low_coverage_regions_length (int): the amount of genes who have a low coverage coding sequence in
        the t2t reference genome. Used for the barplot
        GRCh38_unique_low_coverage_regions_genes_length (int): the amount of genes who have a low coverage coding
        sequence in the GRCh38 reference genomes. Used for the barplot
        common_low_coverage_regions_genes (list): list with the gene_ids who have low coverage coding sequence in both
        reference genomes. Contains duplicates. Used for writing to the file.
        t2t_unique_low_coverage_regions_genes (list): list with the gene_ids who have low coverage coding sequence in
        the T2T reference genome. Contains duplicates. Used for writing to the file.
        GRCh38_unique_low_coverage_regions_genes (list): list with the gene_ids who have low coverage coding sequence in
        the GRCh38 reference genome. Contains duplicates. Used for writing to the file.
        common_low_coverage_regions_genes_dataframe (dataframe): dataframe containing the rows with the gene_ids who
        have low coverage coding sequence in both reference genomes. Contains duplicates. Used for the processing.
        t2t_unique_low_coverage_regions_dataframe (dataframe): dataframe containing the rows with the gene_ids who
        have low coverage coding sequence in the T2T reference genomes. Contains duplicates. Used for the processing.
        GRCh38_unique_low_coverage_regions_dataframe: dataframe containing the rows with the gene_ids who
        have low coverage coding sequence in the GRCh38 reference genomes. Contains duplicates. Used for the processing


    """
    # For each gene_id in the dataframe, checks if it exists in the other dataframe and adds a True/False column to
    # indicate the result
    t2t_dataframe["Gene_id_in_GRCh38"] = t2t_dataframe["Gene_id"].isin(hg38_dataframe["Gene_id"])
    hg38_dataframe["Gene_id_in_T2T"] = hg38_dataframe["Gene_id"].isin(t2t_dataframe["Gene_id"])

    print("Comparing the GeneIDs", "\n")

    # Creates a dataframe with gene_ids where the newly added column is True
    # This dataframe includes gene_ids with coding sequences in low coverage regions for both reference genomes
    common_low_coverage_regions_genes_dataframe = t2t_dataframe[t2t_dataframe["Gene_id_in_GRCh38"]]
    common_low_coverage_regions_genes = list(common_low_coverage_regions_genes_dataframe["Gene_id"])
    common_low_coverage_regions_genes_length = len(list(common_low_coverage_regions_genes_dataframe["Gene_id"]
                                                        .unique()))

    # Creates a dataframe with gene_ids where newly added column in the T2T dataframe is False
    # This dataframe includes gene_ids with coding sequences in low coverage region that are in T2T only
    t2t_unique_low_coverage_regions_dataframe = t2t_dataframe[~t2t_dataframe["Gene_id_in_GRCh38"]]
    t2t_unique_low_coverage_regions_genes = list(t2t_unique_low_coverage_regions_dataframe["Gene_id"])
    t2t_unique_low_coverage_regions_length = len(list(t2t_unique_low_coverage_regions_dataframe["Gene_id"]
                                                      .unique()))


    # Creates a dataframe with gene_ids where newly added column in the GRCh38 dataframe is False
    # This dataframe includes gene_ids with coding sequences in low coverage region that are in GRCh38 only
    GRCh38_unique_low_coverage_regions_dataframe = hg38_dataframe[~hg38_dataframe["Gene_id_in_T2T"]]
    GRCh38_unique_low_coverage_regions_genes = list(GRCh38_unique_low_coverage_regions_dataframe["Gene_id"])
    GRCh38_unique_low_coverage_regions_genes_length = len(list(GRCh38_unique_low_coverage_regions_dataframe["Gene_id"]
                                                               .unique()))

    return (common_low_coverage_regions_genes_length, t2t_unique_low_coverage_regions_length,
            GRCh38_unique_low_coverage_regions_genes_length, common_low_coverage_regions_genes,
            t2t_unique_low_coverage_regions_genes, GRCh38_unique_low_coverage_regions_genes,
            t2t_unique_low_coverage_regions_dataframe, GRCh38_unique_low_coverage_regions_dataframe,
            common_low_coverage_regions_genes_dataframe)


def compare_with_SD(arguments, t2t_unique_low_coverage_regions_dataframe, GRCh38_unique_low_coverage_regions_dataframe,
                    common_low_coverage_regions_genes_dataframe):
    """
    This function identifies overlaps between low coverage regions in the CDS of genes in T2T, GRCh38, or both, and
    the segmental duplications in their respective reference genomes.

    :param:
        arguments.SD_T2T (file): the file containing the segmental duplications in the T2T reference genome.
        arguments.SD_GRCh38 (file): the file containing the segmental duplications in the GRCh38 reference genome.
        t2t_unique_low_coverage_regions_dataframe (dataframe): the dataframe containing the gene_IDs which only have
        low coverage regions in coding sequence in the T2T genome.
        GRCh38_unique_low_coverage_regions_dataframe (dataframe): the dataframe containing the gene_IDs which only have
        low coverage regions in coding sequence in the GRCh38 genome.
        common_low_coverage_regions_genes_dataframe (dataframe): the dataframe containing the gene_IDs which have
        low coverage regions in coding sequence in both genomes.

    :return:
        filtered_t2t_dataframe_length (int): integer showing the amount of gene_ids which have low coverage regions in
        coding sequences and segmental duplications in the T2T reference genome.
        filtered_hg38_dataframe_length (int): integer showing the amount of gene_ids which have low coverage regions in
        coding sequences and segmental duplications in the GRCh38 reference genome.
        unique_common_SD_overlap_gen (int): integer showing the amount of gene_ids which have low coverage regions
        located in coding sequences and segmental duplications in both reference genomes.
        non_overlapping_t2t_genes (dataframe): dataframe containing the gene_IDs whose low coverage coding sequence
        doesn't overlap with a segmental duplication.
        non_overlapping_GRCh38_genes: dataframe containing the gene_IDs whose low coverage coding sequence
        doesn't overlap with a segmental duplication.
    """
    column_names = ["Chromosome", "Start", "End"]

    # Reads the files with the segmental duplications and turns them into dataframes, also adds comments
    t2t_Segmental_Duplications = pd.read_csv(arguments.SD_T2T, sep="\t", encoding="utf-8", names=column_names)
    hg38_Segmental_Duplications = pd.read_csv(arguments.SD_GRCh38, sep="\t", encoding="utf-8", names=column_names)

    print("Finding overlap between the segmental duplications and the low coverage regions", "\n")

    def finding_overlaps_low_coverage_SD(dataframe_1, dataframe_2, suffix1="_df1", suffix2="_df2"):
        """
        This function merges the two different dataframes on the chromosome with the data which overlaps, so no NaN
        values. It then gets the rows where the start of the coding sequence is higher than the start of a
        segmental duplication and lower than the end of a segmental duplication or where the end of the coding sequence
        is higher than the start of a segmental duplication and lower than the end of a segmental duplication.

        :param:
            dataframe_1 (dataframe): the dataframe containing the low coverage regions with their gene_ids.
            dataframe_2 (dataframe): the dataframe containing the segmental duplications.
            suffix1 (string): the suffix of the first dataframe. Gets added to the column names.
            suffix2 (string): the suffix of the second dataframe. Gets added to the column names.

        :return:
            overlap_dataframe (dataframe): the merged dataframe containing the low coverage regions who fully or
            partially overlap with a segmental duplication.
        """
        # Merges the first dataframe with the second on the chromosome with only the data that overlaps
        merged_dataframe = pd.merge(dataframe_1, dataframe_2, on="Chromosome", how="inner", suffixes=(suffix1, suffix2))

        # Gets the rows where the start of the coding sequence is lower or equal to the end of the segmental duplication
        # (if this is not the case then the coding sequence is to the right of the segmental duplication) and where
        # the end of the coding sequence is higher or equal to the start of the segmental duplication (if this is not
        # the case then the coding sequence is to the left of the segmental duplication). These obtained rows overlap
        # with a segmental duplication
        overlap_dataframe = merged_dataframe[
            (merged_dataframe[f"Start{suffix1}"] <= merged_dataframe[f"End{suffix2}"]) &
            (merged_dataframe[f"End{suffix1}"] >= merged_dataframe[f"Start{suffix2}"])
        ]

        return overlap_dataframe

    # Gives the dataframe with the low coverage regions and the dataframe with the segmental duplications in the T2T
    # reference genome to the finding_overlaps_low_coverage_SD function. So, that it can find the low coverage regions
    # which overlap with the segmental duplications and puts them into a dataframe
    filtered_t2t_df = finding_overlaps_low_coverage_SD(t2t_unique_low_coverage_regions_dataframe,
                                                       t2t_Segmental_Duplications,
                                                       suffix1="_low_coverage_region",
                                                       suffix2="_SD")

    # This merges the dataframe with all the genesIDs in GRCh38 with the dataframe containing the geneIDS that overlap
    # with a segmental duplication. Made to get the geneIDs which don't overlap with a segmental duplication
    merged_T2T = pd.merge(t2t_unique_low_coverage_regions_dataframe,
                             filtered_t2t_df[['Chromosome', 'Start_low_coverage_region', 'End_low_coverage_region',
                                               'Gene_id']],
                             left_on=['Chromosome', 'Start', 'End', 'Gene_id'],
                             right_on=['Chromosome', 'Start_low_coverage_region', 'End_low_coverage_region', 'Gene_id'],
                             how='left',
                             indicator=True)

    # Gets the rows of the genes whose coding sequence don't overlap with a segmental duplication by looking which
    # rows are only in the all geneIDs dataframe and dropping which are in both dataframes (all and overlapping with
    # a segmental duplication)
    non_overlapping_t2t_genes = merged_T2T[merged_T2T['_merge'] == 'left_only'].drop('_merge', axis=1)


    # Gives the dataframe with the low coverage regions and the dataframe with the segmental duplications in the GRCh38
    # reference genome to the finding_overlaps_low_coverage_SD function. So, that it can find the low coverage regions
    # which overlap with the segmental duplications and puts them into a dataframe
    filtered_hg38_df = finding_overlaps_low_coverage_SD(GRCh38_unique_low_coverage_regions_dataframe,
                                                        hg38_Segmental_Duplications,
                                                        suffix1="_low_coverage_region",
                                                        suffix2="_SD")

    # This merges the dataframe with all the genesIDs in GRCh38 with the dataframe containing the geneIDS that overlap
    # with a segmental duplication. Made to get the geneIDs which don't overlap with a segmental duplication
    merged_GRCh38 = pd.merge(GRCh38_unique_low_coverage_regions_dataframe,
                      filtered_hg38_df[['Chromosome', 'Start_low_coverage_region', 'End_low_coverage_region',
                                        'Gene_id']],
                      left_on=['Chromosome', 'Start', 'End', 'Gene_id'],
                      right_on=['Chromosome', 'Start_low_coverage_region', 'End_low_coverage_region', 'Gene_id'],
                      how='left',
                      indicator=True)

    # Gets the rows of the genes whose coding sequence don't overlap with a segmental duplication by looking which
    # rows are only in the all geneIDs dataframe and dropping which are in both dataframes (all and overlapping with
    # a segmental duplication)
    non_overlapping_GRCh38_genes = merged_GRCh38[merged_GRCh38['_merge'] == 'left_only'].drop('_merge', axis=1)

    # Gives the dataframe with the low coverage regions and Gene_ids in both reference genomes and
    # the dataframe with the segmental duplications in the T2T reference genome to the
    # finding_overlaps_low_coverage_SD function. So, that it can find the low coverage regions which overlap with
    # the segmental duplications and puts them into a dataframe
    merged_dataframe_common_t2t = finding_overlaps_low_coverage_SD(common_low_coverage_regions_genes_dataframe,
                                                                   t2t_Segmental_Duplications,
                                                                   suffix1="_common_low_coverage",
                                                                   suffix2="_SD_T2T")

    # Changes column names to prevent overwriting the previous gotten low coverage regions
    merged_dataframe_common_t2t.columns = ["Chromosome", "Start", "End", "Gene_id", "Gene_id_in_GRCh38", "Start_T2T",
                                           "End_T2T"]

    # Gives the merged dataframe with the low coverage regions and Gene_ids which overlap with Segmental duplications
    # from the T2T reference genome and the dataframe with the segmental duplications in the GRCh38 reference genome to
    # the finding_overlaps_low_coverage_SD function. So, that it can find the low coverage regions which overlap with
    # the segmental duplications and puts them into a dataframe
    merged_dataframe_with_Hg38 = finding_overlaps_low_coverage_SD(merged_dataframe_common_t2t,
                                                                  hg38_Segmental_Duplications,
                                                                  suffix1="_common_low_coverage",
                                                                  suffix2="_SD_GRCh38")

    # Combines the unique SD overlapping gene_IDS from both merged dataframes without a gene_id getting duplicated. It
    # counts a geneID overlapping if it overlaps with a segmental duplication in one or both reference genomes.
    common_low_coverage_SD_overlap = numpy.unique(numpy.concatenate([merged_dataframe_common_t2t["Gene_id"]
                                                                    .unique(),
                                                                    merged_dataframe_with_Hg38["Gene_id"].unique()]
                                                                    ,axis=0))

    # Gets the amount of unique gene_IDS which have low coverage regions in coding sequences which overlap with a
    # segmental duplication
    filtered_t2t_dataframe_length = len(filtered_t2t_df["Gene_id"].unique())
    filtered_hg38_dataframe_length = len(filtered_hg38_df["Gene_id"].unique())
    unique_common_SD_overlap_genes = len(common_low_coverage_SD_overlap)

    return (filtered_t2t_dataframe_length, filtered_hg38_dataframe_length, unique_common_SD_overlap_genes,
            non_overlapping_t2t_genes, non_overlapping_GRCh38_genes)



def make_barplot(t2t_unique_low_coverage_regions_length,
                 GRCh38_unique_low_coverage_regions_length, filtered_t2t_dataframe_length,
                 filtered_hg38_dataframe_length):
    """
    This function makes the barplot where the length of the lists made in the previous function are plotted against
    each other.

    :param:
        t2t_unique_low_coverage_regions_length (int): integer showing gene_ids with coding sequences located
        exclusively in low coverage regions of T2T.
        GRCh38_unique_low_coverage_regions_length (int): integer showing gene_ids with coding sequences located
        exclusively in low coverage regions of GRCh38.
        filtered_t2t_dataframe_length (int): integer showing the amount of gene_ids which have low coverage regions in
        coding sequences and segmental duplications in the T2T reference genome.
        filtered_hg38_dataframe_length (int): integer showing the amount of gene_ids which have low coverage regions in
        coding sequences and segmental duplications in the GRCh38 reference genome.

    :return:
        The barplot is saved as Coding_Sequences_in_low_coverage_regions.png
    """
    print("Making the barplot", "\n")

    # Makes matplotlib work on the server.
    matplotlib.use("Agg")

    # Makes subplots
    fig, ax = plt.subplots()

    # The labels for the x and the legend
    xlables = ["T2T: lower coverage \n Hg38: higher coverage ", "T2T: higher coverage \n Hg38: lower coverage "]


    # Makes the barplot that shows the amount of genes which have a coding sequence overlap with a low coverage region
    bars = ax.bar(xlables,
                  [t2t_unique_low_coverage_regions_length,
                   GRCh38_unique_low_coverage_regions_length],
                  alpha=0.5,
                  label="no segmental duplication overlap",
                  )

    # Makes an additional barplot containing the amount of genes whose low coverage coding sequence overlap with a
    # segmental duplication
    SD_bars = ax.bar(xlables, [filtered_t2t_dataframe_length,
                     filtered_hg38_dataframe_length], color=["black", "black"], alpha=0.9,
                     label="segmental duplication overlap")

    # Adds amount of total genes per category to their respective bar plots by calculating the midpoint of each bar and
    # placing the values there.
    for bar in bars:
        height = bar.get_height()
        ax.annotate("{}".format(height),
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0,5),
                    textcoords="offset points",
                    ha="center",
                    va="center",
                    fontsize="medium")

    # Adds amount of genes whose coding sequence overlap with a segmental duplication to their respective bar plots by
    # calculating the midpoint of each bar and placing the values there.
    for bar in SD_bars:
        height = bar.get_height()
        ax.annotate("{}".format(height),
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 5),
                    textcoords="offset points",
                    ha="center",
                    va="center")

    # Makes and sets the y_label
    ax.set_ylabel("Amount of genes")

    # Makes and sets the title
    plt.title("Amount of genes with coding sequences in low coverage regions")
    plt.legend(fontsize="small")

    # Saves the barplot into a png file
    plt.savefig("Coding_Sequences_in_low_coverage_regions.png")
    print("The barplot is made and saved as Coding_Sequences_in_low_coverage_regions.png")
    # plt.show()


def write_to_file(arguments, common_genes, t2t_unique_genes, GRCh38_unique_genes, non_overlapping_t2t_genes,
                  non_overlapping_GRCh38_genes):
    """
    This function writes the genes of each dataframe to their respective file.

    :param:
        arguments.Low_coverage_genes_shared (txt file): empty text file where the geneIDs of the genes which have CDS
        in a low coverage region in both reference genome will be written to.
        arguments.Low_coverage_genes_T2T (txt file): empty text file where the geneIDs of the genes which have CDS
        in a low coverage region in T2T wil lbe written to.
        arguments.Low_coverage_genes_GRCh38 (txt file): empty text file where the geneIDs of the genes which have CDS
        in a low coverage region in GRCh38 will be written to.
        common_genes (list): list containing the geneIDs of the genes which have CDS
        in a low coverage region in both reference genomes.
        t2t_unique_genes (list): list containing the geneIDs of the genes which have a CDS
        in a low coverage region in T2T.
        GRCh38_unique_genes (list): list containing the geneIDs of the genes which have a CDS
        in a low coverage region in GRCh38.
        non_overlapping_t2t_genes (dataframe): dataframe containing the genes whose low coverage CDS don't overlap
        with a segmental duplication
        non_overlapping_GRCh38_genes (dataframe): dataframe containing the genes whose low coverage CDS don't overlap
        with a segmental duplication

    :return:
        nothing
    """
    # Opens the files and writes the list to that specific file
    with open(arguments.Low_coverage_genes_shared, "w") as shared_file:
        for gene in common_genes:
            shared_file.write("{}\n".format(gene))

    with open(arguments.Low_coverage_genes_T2T, "w") as t2t_file:
        for gene in t2t_unique_genes:
            t2t_file.write("{}\n".format(gene))

    with open(arguments.Low_coverage_genes_GRCh38, "w") as GRCh38_file:
        for gene in GRCh38_unique_genes:
            GRCh38_file.write("{}\n".format(gene))

    # Opens the empty bed file and writes the dataframes to that file
    non_overlapping_t2t_genes.to_csv(arguments.non_SD_overlapping_genes_T2T, index=False, sep="\t")

    non_overlapping_GRCh38_genes.to_csv(arguments.non_SD_overlapping_genes_Hg38, index=False, sep="\t")

def main(args):
    t2t_dataframe, hg38_dataframe = get_files(args)

    (common_genes_length, t2t_unique_length, GRCh38_unique_length, common_genes, t2t_unique_genes,
     GRCh38_unique_genes, t2t_unique_dataframe, GRCh38_unique_dataframe, common_dataframe) = sort_files(t2t_dataframe,
                                                                                                        hg38_dataframe)

    (filtered_t2t_dataframe_length, filtered_hg38_dataframe_length, unique_common_SD_overlap_genes,
     non_overlapping_t2t_genes, non_overlapping_GRCh38_genes) = compare_with_SD(args,
                                                                                t2t_unique_dataframe,
                                                                                GRCh38_unique_dataframe,
                                                                                common_dataframe)

    make_barplot(
                 t2t_unique_length,
                 GRCh38_unique_length,
                 filtered_t2t_dataframe_length,
                 filtered_hg38_dataframe_length,
                 )

    write_to_file(args,
                  common_genes,
                  t2t_unique_genes,
                  GRCh38_unique_genes,
                  non_overlapping_t2t_genes,
                  non_overlapping_GRCh38_genes)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("T2T_intersected_CDS_Only",
                        help="Path to the Intersected bed file from the T2T RefSeq and categorical files containing the"
                             " genes whose coding sequences are in low coverage regions",
                        metavar="Intersected bed file from the T2T RefSeq and categorical files")
    parser.add_argument("GRCh38_intersected_CDS_only",
                        help="Intersected bed file from the GRCh38 RefSeq and categorical files containing the genes "
                                "whose coding sequences are in low coverage regions",
                        metavar="Intersected bed file")
    parser.add_argument("Low_coverage_genes_shared",
                        help="Path to a empty txt file where the gene_ids with low coverage coding sequences in both "
                             "reference genomes will be sent to")
    parser.add_argument("Low_coverage_genes_T2T",
                        help="Path to a empty txt file where the gene_ids with low coverage coding sequences in the "
                        "T2T reference genomes will be sent to")
    parser.add_argument("Low_coverage_genes_GRCh38",
                        help="Path to a empty txt file where the gene_ids with low coverage coding sequences in the "
                             "GRCh38 reference genomes will be sent to")
    parser.add_argument("SD_T2T",
                        help="Path to the file containing the known segmental duplications in the T2T reference genome")
    parser.add_argument("SD_GRCh38",
                        help="Path to the file containing the segmental duplications in the GRCh38 reference genome")
    parser.add_argument("non_SD_overlapping_genes_T2T",
                        help="An empty file where the low coverage genes which don't overlap with an "
                             "segmental duplication in the T2T reference genome will be sent to")
    parser.add_argument("non_SD_overlapping_genes_Hg38",
                        help="An empty file where the low coverage genes which don't overlap with an "
                             "segmental duplication in the GRCh38 reference genome will be sent to")

    args = parser.parse_args()
    main(args)