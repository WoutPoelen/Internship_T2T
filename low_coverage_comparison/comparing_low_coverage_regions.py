import argparse
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd


def generate_dataframes(args):


    t2t_dataframe = pd.read_csv(args.T2T_BED_file, sep="\t", encoding="utf-8")
    hg38_dataframe = pd.read_csv(args.GRCh38_BED_file, sep="\t", encoding="utf-8")
    liftover_dataframe = pd.read_csv(args.Liftover_BED_file, sep="\t", encoding="utf-8")

    print("Reading the files")

    t2t_dataframe.columns = ["Chromosome", "Start", "End", "mean_coverage"]
    hg38_dataframe.columns = ["Chromosome", "Start", "End", "mean_coverage"]
    liftover_dataframe.columns = ["Chromosome", "Start", "End", "mean_coverage", "succesful_liftover"]

    print("Processing the files")
    hg38_dataframe = hg38_dataframe[(hg38_dataframe["End"] <= 5000000) & (hg38_dataframe["Chromosome"] == "chr1") &
                                    (hg38_dataframe["mean_coverage"] <= 10)]
    t2t_dataframe = t2t_dataframe[(t2t_dataframe["End"] <= 5000000) & (t2t_dataframe["Chromosome"] == "chr1") &
                                  (t2t_dataframe["mean_coverage"] <= 10)]

    return t2t_dataframe, hg38_dataframe, liftover_dataframe


def plot_low_coverage(dataframe_t2t, dataframe_hg38, dataframe_liftover):
    # Makes the plotting work on the server
    matplotlib.use("Agg")

    print("Generating the plot")

    dataframe_t2t["Reference genome"] = "T2T"
    dataframe_hg38["Reference genome"] = "GRCh38"
    dataframe_liftover["Reference genome"] = "liftover_T2T"

    figs, ax = plt.subplots(figsize=(10, 5))
    plt.scatter(dataframe_t2t["Start"], dataframe_t2t["Reference genome"], s=10)
    plt.scatter(dataframe_liftover["Start"], dataframe_liftover["Reference genome"], s=10)
    plt.scatter(dataframe_hg38["Start"], dataframe_hg38["Reference genome"], s=10)

    plt.title("P3-D10 chr1")

    plt.xlabel("Position (Mbp)")

    plt.savefig("low_coverage_comparison.png")

    print("The low coverage plot has been successfully generated and save as low_coverage_comparison.png")
    # plt.show()


def main(args):
    t2t_dataframe, hg38_dataframe, liftover_dataframe = generate_dataframes(args)
    plot_low_coverage(t2t_dataframe, hg38_dataframe, liftover_dataframe)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("T2T_BED_file",
                        help="Path to the filtered BED file from process_regions_file.py containing "
                             "the regions of 500 bp and the average coverage lower than 10 from the T2T file")
    parser.add_argument("GRCh38_BED_file",
                        help="Path to the filtered BED file from process_regions_file.py containing the regions "
                             "of 500 bp and the average coverage lower than 10 from the GRCh38 file")
    parser.add_argument("Liftover_BED_file",
                        help="Path to the from T2T lifted over BED file containing the coordinates for the locations"
                             "of the lifted over regions on the GRCh38 reference genome")
    args = parser.parse_args()
    main(args)