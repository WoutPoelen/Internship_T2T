library(karyoploteR)
library(GenomicRanges)

data_GRCh38 <- toGRanges(data.frame(chr=GRCh38_entire_genome_P3.D10_low_coverage["V1"], 
                                  start=GRCh38_entire_genome_P3.D10_low_coverage["V2"], 
                                  end=GRCh38_entire_genome_P3.D10_low_coverage["V3"],
                                  y=GRCh38_entire_genome_P3.D10_low_coverage["V4"]))

plot_params <- getDefaultPlotParams(plot.type = 1)

plot_params$data1height <-  200
plot_params$topmargin <- 250
plot_params$leftmargin <- 0.06
plot_params$ideogramlateralmargin <- 0.05
plot_params$data.panel.margin <- 0.1

kp_GRCh38 <- plotKaryotype(genome = "hg38", plot.params = plot_params, 
                           plot.type = 1)
kpAddBaseNumbers(kp_GRCh38)
kpPlotDensity(kp_GRCh38, data=data_GRCh38, data.panel = 1, window.size = 5000,
              col="black")
kpAddMainTitle(kp_GRCh38, "Low coverage regions GRCh38", col="black")


