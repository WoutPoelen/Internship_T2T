import argparse
import matplotlib
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
    # Turns the two bed files into pandas dataframes
    t2t_dataframe = pd.read_csv(arguments.T2T_intersected_CDS_Only, sep="\t", encoding="utf-8")
    hg38_dataframe = pd.read_csv(arguments.GRCh38_intersected_CDS_only, sep="\t", encoding="utf-8")

    # Adds columns to the dataframes
    hg38_dataframe.columns = ["Chromosome", "Start", "End", "Gene_id"]
    t2t_dataframe.columns = ["Chromosome", "Start", "End", "Gene_id"]

    # Ensure each gene_id with overlapping CDS appears only once in the dataframe to avoid duplicate counts.
    t2t_dataframe = t2t_dataframe.drop_duplicates(subset="Gene_id")
    hg38_dataframe = hg38_dataframe.drop_duplicates(subset="Gene_id")

    # Add _Y to the end of the genes on the Y chromosome, because there are genes that are on both X and Y chromosome,
    # but only on the Y chromosome is it low coverage
    t2t_condition = (t2t_dataframe["Chromosome"] == "chrY")
    hg38_condition = (hg38_dataframe["Chromosome"] == "chrY")
    t2t_dataframe.loc[t2t_condition, "Gene_id"] = t2t_dataframe.loc[t2t_condition, "Gene_id"].astype(str) + "_Y"
    hg38_dataframe.loc[hg38_condition, "Gene_id"] = hg38_dataframe.loc[hg38_condition, "Gene_id"].astype(str) + "_Y"

    return t2t_dataframe, hg38_dataframe

def sort_files(t2t_dataframe, hg38_dataframe):
    """
    This function takes the dataframes made in the previous function and compares the gene_ids between the dataframes to
    see which gene_ids have a CDS which overlaps with a low coverage region.

    Adds the gene_id based on the comparison result:
        - If the gene_id is in both dataframes, it is added to common_genes_length.
        - If the gene_id is only in the T2T dataframe, it is added to t2t_unique_length.
        - If the gene_id is only in the GRCh38 dataframe, it is added to GRCh38_unique_length.

    The lengths of common_genes_length, t2t_unique_length, and GRCh38_unique_length are then calculated and saved.

    :param:
        t2t_dataframe (dataframe): Dataframe containing the chromosome, start, end and Gene_id of the CDS regions in
        T2T with every Gene_id only appearing once.
        hg38_dataframe: Dataframe containing the chromosome, start, end and Gene_id of the CDS regions in
        GRCh38 with every Gene_id only appearing once.

    :return:
        common_genes_length (dataframe): dataframe containing the gene_ids which have coding sequences in low coverage regions
        in both T2T and GRCh38.
        t2t_unique_length (dataframe): dataframe containing gene_ids with coding sequences located exclusively in low coverage
        regions of T2T.
        GRCh38_unique_length (dataframe):dataframe containing gene_ids with coding sequences located exclusively in low
        coverage regions of GRCh38.
    """
    # For each gene_id in the dataframe, checks if it exists in the other dataframe and adds a True/False column to
    # indicate the result
    t2t_dataframe["Gene_id_in_GRCh38"] = t2t_dataframe["Gene_id"].isin(hg38_dataframe["Gene_id"])
    hg38_dataframe["Gene_id_in_T2T"] = hg38_dataframe["Gene_id"].isin(t2t_dataframe["Gene_id"])

    # Creates a dataframe with gene_ids where the newly added column is True
    # This dataframe includes gene_ids with coding sequences in low coverage regions for both reference genomes
    common_genes_dataframe = t2t_dataframe[t2t_dataframe["Gene_id_in_GRCh38"]]
    common_genes = list(common_genes_dataframe["Gene_id"])
    common_genes_length = len(list(common_genes_dataframe["Gene_id"]))

    # Creates a dataframe with gene_ids where newly added column in the T2T dataframe is False
    # This dataframe includes gene_ids with coding sequences in low coverage region that are in T2T only
    t2t_unique_dataframe = t2t_dataframe[~t2t_dataframe["Gene_id_in_GRCh38"]]
    t2t_unique_genes = list(t2t_unique_dataframe["Gene_id"])
    t2t_unique_length = len(list(t2t_unique_dataframe["Gene_id"]))

    # Creates a dataframe with gene_ids where newly added column in the GRCh38 dataframe is False
    # This dataframe includes gene_ids with coding sequences in low coverage region that are in GRCh38 only
    GRCh38_unique_dataframe = hg38_dataframe[~hg38_dataframe["Gene_id_in_T2T"]]
    GRCh38_unique_genes = list(GRCh38_unique_dataframe["Gene_id"])
    GRCh38_unique_length = len(list(GRCh38_unique_dataframe["Gene_id"]))

    return (common_genes_length, t2t_unique_length, GRCh38_unique_length, common_genes, t2t_unique_genes,
            GRCh38_unique_genes)


def make_barplot(common_genes_length, t2t_unique_length, GRCh38_unique_length):
    """
    This function makes the barplot where the length of the lists made in the previous function are plotted against
    each other.

    :param:
        common_genes_length (int): integer showing the amount ofe gene_ids which have coding sequences
        in low coverage regions in both T2T and GRCh38.
        t2t_unique (int): integer showing gene_ids with coding sequences located exclusively in low coverage
        regions of T2T.
        GRCh38_unique (int): integer showing gene_ids with coding sequences located exclusively in low
        coverage regions of GRCh38.

    :return:
        The barplot is saved as Coding_Sequences_in_low_coverage_regions.png
    """
    # Makes matplotlib work on the server.
    matplotlib.use("Agg")

    # Makes subplots
    fig, ax = plt.subplots()

    # Makes the barplot
    bars = ax.bar(["shared", "T2T only", "GRCh38 only"], [common_genes_length, t2t_unique_length,
                                                          GRCh38_unique_length], color=["blue", "green", "red"])

    # Adds variable values to their respective bar plots by calculating the midpoint of each bar and placing the values
    # there.
    for bar in bars:
        height = bar.get_height()
        ax.annotate("{}".format(height),
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0,5),
                    textcoords="offset points",
                    ha="center",
                    va="center")

    # Makes and sets the y_label
    ax.set_ylabel("Amount of genes")

    # Makes and sets the title
    plt.title("Amount of genes with coding sequences in low coverage regions")

    # Saves the barplot into a png file
    plt.savefig("Coding_Sequences_in_low_coverage_regions.png")
    # plt.show()


def write_to_file(arguments, common_genes, t2t_unique_genes, GRCh38_unique_genes):
    """

    :param:
        arguments.Low_coverage_genes_shared (txt file): empty text file where the geneIDs of the genes which have CDS
        in a low coverage region in both reference genome will be written to.
        arguments.Low_coverage_genes_T2T (txt file): empty text file where the geneIDs of the genes which have CDS
        in a low coverage region in T2T wil lbe written to.
        arguments.Low_coverage_genes_GRCh38 (txt file): empty text file where the geneIDs of the genes which have CDS
        in a low coverage region in GRCh38 will be written to.
        common_genes (list): list containing the geneIDs of the genes which have CDS
        in a low coverage region in both reference genomes.
        t2t_unique_genes (list): list containing the geneIDs of the genes which have coding CDS
        in a low coverage region in T2T.
        GRCh38_unique_genes (list): list containing the geneIDs of the genes which have coding CDS
        in a low coverage region in GRCh38.
    :return:
    """
    with open(arguments.Low_coverage_genes_shared, "w") as shared_file:
        for gene in common_genes:
            shared_file.write("{}\n".format(gene))

    with open(arguments.Low_coverage_genes_T2T, "w") as t2t_file:
        for gene in t2t_unique_genes:
            t2t_file.write("{}\n".format(gene))

    with open(arguments.Low_coverage_genes_GRCh38, "w") as GRCh38_file:
        for gene in GRCh38_unique_genes:
            GRCh38_file.write("{}\n".format(gene))



def main(args):
    t2t_dataframe, hg38_dataframe = get_files(args)
    (common_genes_length, t2t_unique_length, GRCh38_unique_length, common_genes, t2t_unique_genes,
     GRCh38_unique_genes) = sort_files(t2t_dataframe, hg38_dataframe)
    make_barplot(common_genes_length, t2t_unique_length, GRCh38_unique_length)
    write_to_file(args, common_genes, t2t_unique_genes, GRCh38_unique_genes)

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

    args = parser.parse_args()
    main(args)