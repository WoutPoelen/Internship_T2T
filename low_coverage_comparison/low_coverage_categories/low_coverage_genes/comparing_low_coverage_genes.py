import argparse
import matplotlib
import pandas as pd
from matplotlib import pyplot as plt


def get_files(arguments):
    # Unzips with gzip and turns the files given by the user into dataframes
    t2t_dataframe = pd.read_csv(arguments.T2T_filtered_categorical, sep="\t", encoding="utf-8")
    hg38_dataframe = pd.read_csv(arguments.GRCh38_filtered_categorical, sep="\t", encoding="utf-8")


    hg38_dataframe.columns = ["Chromosome", "Start", "End", "Gene_id"]
    t2t_dataframe.columns = ["Chromosome", "Start", "End", "Gene_id"]
    
    t2t_dataframe = t2t_dataframe.drop_duplicates(subset="Gene_id")
    hg38_dataframe = hg38_dataframe.drop_duplicates(subset="Gene_id")
    print(t2t_dataframe)
    print(hg38_dataframe)


    # Add _Y to the end of the genes on the Y chromosome, because there are genes that are on both X and Y chromosome,
    # but only on the Y chromosome is it low coverage
    t2t_condition = (t2t_dataframe["Chromosome"] == "chrY")
    hg38_condition = (hg38_dataframe["Chromosome"] == "chrY")
    t2t_dataframe.loc[t2t_condition, "Gene_id"] = t2t_dataframe.loc[t2t_condition, "Gene_id"].astype(str) + "_Y"
    hg38_dataframe.loc[hg38_condition, "Gene_id"] = hg38_dataframe.loc[hg38_condition, "Gene_id"].astype(str) + "_Y"

    return t2t_dataframe, hg38_dataframe

def sort_files(t2t_dataframe, hg38_dataframe):
    # merge_dataframe= pd.merge(t2t_dataframe, hg38_dataframe, on= ["Chromosome", "Start", "End"], how="outer",
    #                           indicator=True)

    t2t_dataframe["Gene_id_in_GRCh38"] = t2t_dataframe["Gene_id"].isin(hg38_dataframe["Gene_id"])
    hg38_dataframe["Gene_id_in_T2T"] = hg38_dataframe["Gene_id"].isin(t2t_dataframe["Gene_id"])

    common_genes_dataframe = t2t_dataframe[t2t_dataframe["Gene_id_in_GRCh38"]]
    common_genes = len(list(common_genes_dataframe["Gene_id"]))
    t2t_unique_dataframe = t2t_dataframe[~t2t_dataframe["Gene_id_in_GRCh38"]]
    t2t_unique = len(list(t2t_unique_dataframe["Gene_id"]))
    GRCh38_unique_dataframe = hg38_dataframe[~hg38_dataframe["Gene_id_in_T2T"]]
    GRCh38_unique = len(list(GRCh38_unique_dataframe["Gene_id"]))

    print("t2t unique: ", t2t_unique, "\n",
          "GRCh38 unique: ", GRCh38_unique, "\n",
          "common genes: ", common_genes)

    return common_genes, t2t_unique, GRCh38_unique


def make_barplot(common_genes, t2t_unique, GRCh38_unique):

    # matplotlib.use("Agg")
    fig, ax = plt.subplots()

    bars = ax.bar(["shared", "T2T unique", "GRCh38 unique"], [common_genes, t2t_unique, GRCh38_unique],
                  color=["blue", "green", "red"])

    for bar in bars:
        height = bar.get_height()
        ax.annotate("{}".format(height),
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0,5),
                    textcoords="offset points",
                    ha="center",
                    va="center")
    ax.set_ylabel("Amount of genes")
    plt.title("Amount of genes in low coverage regions")

    plt.show()

def main(args):
    t2t_dataframe, hg38_dataframe = get_files(args)
    common_genes, t2t_unique, GRCh38_unique = sort_files(t2t_dataframe, hg38_dataframe)
    make_barplot(common_genes, t2t_unique, GRCh38_unique)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("T2T_filtered_categorical",
                        help="Path to the file containing the genes that are in low coverage regions in both"
                             "T2T and GRCh38")
    parser.add_argument("GRCh38_filtered_categorical",
                        help="Path to the file containing the genes that are in low coverage regions in only the"
                             "T2T")

    args = parser.parse_args()
    main(args)