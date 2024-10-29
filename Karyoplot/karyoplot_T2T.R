library(karyoploteR)
library(GenomicRanges)

data_T2T <- toGRanges(data.frame(chr=T2T_low_coverage_HG002["V1"], 
                             start=T2T_low_coverage_HG002["V2"], 
                             end=T2T_low_coverage_HG002["V3"],
                             y=T2T_low_coverage_HG002["V4"]))

custom_genome <- toGRanges(data.frame(chr=chromosome_report_t2t$UCSC.style.name,
                                    start=0, 
                                    end=chromosome_report_t2t$Seq.length))

custom_cytobands <- toGRanges(data.frame(chr=filtered_centromeres_CenSat$V1,
                                         start=filtered_centromeres_CenSat$V2,
                                         end=filtered_centromeres_CenSat$V3))

plot_params <- getDefaultPlotParams(plot.type = 1)

plot_params$data1height <-  200
plot_params$topmargin <- 250
plot_params$leftmargin <- 0.06
plot_params$ideogramlateralmargin <- 0.05
plot_params$data.panel.margin <- 0.1

kp_t2t <- plotKaryotype(genome = custom_genome,plot.type = 1, 
                        chromosomes = "all", plot.params = plot_params)
kpAddBaseNumbers(kp_t2t)
kpPlotRegions(kp_t2t, data = custom_cytobands, data.panel = "ideogram", 
              col = "red")
kpPlotDensity(kp_t2t, data=data_T2T, data.panel = 1, window.size = 5000, 
              col="black")
kpAddMainTitle(kp_t2t, "Low coverage regions T2T-CHM13", col="black")

