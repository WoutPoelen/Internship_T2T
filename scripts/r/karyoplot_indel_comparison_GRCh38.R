library(karyoploteR)
library(GenomicRanges)


GRCh38_indel <- read.table(
  "C:/Users/Z659216/Desktop/BAM_files/svs_files/filtered_centromeres_intersected_GRCh38.bed", 
  header=FALSE, 
  sep="\t",
  stringsAsFactors = FALSE)

total_deletions_GRCh38 <- GRCh38_indel[(GRCh38_indel$V5 == "DEL"), ]
total_insertions_GRCh38 <- GRCh38_indel[(GRCh38_indel$V5 == "INS"), ]

total_deletions_GRCh38_data <- GRCh38_insertion_data <- toGRanges(
  data.frame(chr=total_deletions_GRCh38["V1"], 
             start=total_deletions_GRCh38["V2"], 
             end=total_deletions_GRCh38["V3"],
             y=total_deletions_GRCh38["V4"]))

total_insertions_GRCh38_data <- toGRanges(
  data.frame(chr=total_insertions_GRCh38["V1"], 
             start=total_insertions_GRCh38["V2"], 
             end=total_insertions_GRCh38["V3"],
             y=total_insertions_GRCh38["V4"]))

plot_params <- getDefaultPlotParams(plot.type = 2)

plot_params$data1height <-  100
plot_params$data2height <- 100
plot_params$topmargin <- 250
plot_params$leftmargin <- 0.06
plot_params$ideogramlateralmargin <- 0.05
plot_params$data.panel.margin <- 0.1

kp_GRCh38_indels <- plotKaryotype(genome = "hg38", plot.params = plot_params, 
                                  plot.type = 2)

kpDataBackground(kp_GRCh38_indels, data.panel = 1, col="pink")
kpDataBackground(kp_GRCh38_indels, data.panel = 2, col="lightblue")

kpPlotRegions(kp_GRCh38_indels, data=total_insertions_GRCh38, data.panel = 1,
              col="black")

kpPlotRegions(kp_GRCh38_indels, data=total_deletions_GRCh38, data.panel = 2,
              col="red")

kpAddMainTitle(kp_GRCh38_indels, "Indels in GRCh38")
legend(x="bottomright", fill= c("black", "red"), legend = c("insertions", "deletions"))




