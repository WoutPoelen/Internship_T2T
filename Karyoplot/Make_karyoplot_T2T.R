library(karyoploteR)
library(GenomicRanges)

args <- commandArgs(trailingOnly = TRUE)

file_path_t2t_low_coverage <- args[1]
chromosome_locations <- args[2]
centromeres_locations <- args[3]
path_to_plot <- args[4]

png(path_to_plot)

bed_data <- read.table(file_path_t2t_low_coverage, header = FALSE, sep = "\t")

chromosomes <- read.table(chromosome_locations, header = TRUE, sep = "\t")

centromeres <- read.table(centromeres_locations, header = FALSE, sep = "\t")

colnames(bed_data) <- c("chromosome", "start", "end", "mean_coverage")
colnames(centromeres) <- c("chromosome", "start", "end") 

data_T2T <- toGRanges(data.frame(chr = bed_data$chromosome, 
                                 start = bed_data$start, 
                                 end = bed_data$end,
                                 y = bed_data$mean_coverage))

custom_genome <- toGRanges(data.frame(chr = chromosome_report_t2t$UCSC.style.name,
                                      start = 0, 
                                      end = chromosome_report_t2t$Seq.length))

custom_cytobands <- toGRanges(data.frame(chr = centromeres$start,
                                         start = centromeres$start,
                                         end = centromeres$end))

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

dev.off()