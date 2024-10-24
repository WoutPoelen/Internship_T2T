library(karyoploteR)


args <- commandArgs(trailingOnly = TRUE)

file_path_GRCh38_low_coverage <- args[1]
path_to_plot <- args[2]

png(path_to_plot)

bed_data <- read.table(file_path_GRCh38_low_coverage, header = FALSE, sep = "\t")

colnames(bed_data) <- c("chromosome", "start", "end", "mean_coverage")

data_GRCh38 <- toGRanges(data.frame(chr = bed_data$chromosome, 
                                 start = bed_data$start, 
                                 end = bed_data$end,
                                 y = bed_data$mean_coverage))

plot_params <- getDefaultPlotParams(plot.type = 1)

plot_params$data1height <-  200
plot_params$topmargin <- 250
plot_params$leftmargin <- 0.06
plot_params$ideogramlateralmargin <- 0.05
plot_params$data.panel.margin <- 0.1

kp_GRCh38 <- plotKaryotype(genome = "hg38",plot.type = 1, 
                        chromosomes = "all", plot.params = plot_params)
kpAddBaseNumbers(kp_GRCh38)
kpPlotDensity(kp_GRCh38, data=data_GRCh38, data.panel = 1, window.size = 5000, 
              col="black")
kpAddMainTitle(kp_GRCh38, "Low coverage regions GRCh38", col="black")

def.off()

