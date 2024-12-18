library(karyoploteR)
library(GenomicRanges)

T2T_indel <- read.table(
  #"C:/Users/Z659216/Desktop/BAM_files/svs_files/filtered_centromere_intersected_T2T_indel.bed", 
  "C:/Users/Z659216/Desktop/BAM_files/multiple_sample/P1-B3.bed",
  header=FALSE, 
  sep="\t",
  stringsAsFactors = FALSE)

total_deletions_T2T <- T2T_indel[(T2T_indel$V5 == "DEL"), ]
total_insertions_T2T <- T2T_indel[(T2T_indel$V5 == "INS"), ]

total_T2T_insertion_data <- toGRanges(data.frame(chr=total_insertions_T2T["V1"], 
                                           start=total_insertions_T2T["V2"], 
                                           end=total_insertions_T2T["V3"],
                                           y=total_insertions_T2T["V4"]))

total_T2T_deletion_data <- toGRanges(data.frame(chr=total_deletions_T2T["V1"], 
                                          start=total_deletions_T2T["V2"], 
                                          end=total_deletions_T2T["V3"],
                                          y=total_deletions_T2T["V4"]))

custom_genome <- toGRanges(data.frame(chr=chromosome_report_t2t$UCSC.style.name,
                                      start=0, 
                                      end=chromosome_report_t2t$Seq.length))

custom_cytobands <- toGRanges(data.frame(chr=filtered_centromeres_CenSat$V1,
                                         start=filtered_centromeres_CenSat$V2,
                                         end=filtered_centromeres_CenSat$V3))

plot_params <- getDefaultPlotParams(plot.type = 2)

plot_params$data1height <-  100
plot_params$data2height <- 100
plot_params$topmargin <- 250
plot_params$leftmargin <- 0.06
plot_params$ideogramlateralmargin <- 0.05
plot_params$data.panel.margin <- 0.1

kp_t2t_indel <- plotKaryotype(genome = custom_genome,plot.type = 2, 
                              chromosomes = "all", plot.params = plot_params)

kpDataBackground(kp_t2t_indel, data.panel = 1, col="pink")
kpDataBackground(kp_t2t_indel, data.panel = 2, col="lightblue")

kpPlotRegions(kp_t2t_indel, data = custom_cytobands, data.panel = "ideogram", 
              col = "blue")
kpPlotRegions(kp_t2t_indel, data=total_T2T_insertion_data, data.panel = 1, 
              col="black")
kpPlotRegions(kp_t2t_indel, data=total_T2T_deletion_data, data.panel = 2, 
              col="red")
kpAddMainTitle(kp_t2t_indel, "Indels in T2T", col="black")

legend(x="bottomright", fill= c("black", "red"), legend = c("insertions", "deletions"))
