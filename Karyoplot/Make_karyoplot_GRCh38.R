library(karyoploteR)
library(GenomicRanges)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

file_path_GRCh38_low_coverage <- args[1]
path_to_plot <- args[2]

pdf(path_to_plot, width = 500, height = 100)

bed_data <- read.table(file_path_GRCh38_low_coverage, header = FALSE, sep = "\t")


colnames(bed_data) <- c("chromosome", "start", "end", "mean_coverage")

normal_chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

filtered_data <- bed_data %>% filter(bed_data$chromosome %in% normal_chromosomes)

data_GRCh38 <- toGRanges(data.frame(chr = normal_chromosomes, 
                                 start = filtered_data$start, 
                                 end = filtered_data$end,
                                 y = filtered_data$mean_coverage))

plot_params <- getDefaultPlotParams(plot.type = 1)

plot_params$data1height <-  200
plot_params$topmargin <- 250
plot_params$leftmargin <- 0.06
plot_params$ideogramlateralmargin <- 0.05
plot_params$data.panel.margin <- 0.1

kp_GRCh38 <- plotKaryotype(genome = "hg38", plot.type = 1, 
                        chromosomes = "all") #, plot.params = plot_params)
kpAddBaseNumbers(kp_GRCh38)
kpPlotDensity(kp_GRCh38, data=data_GRCh38, data.panel = 1, window.size = 5000, 
              col="black")
kpAddMainTitle(kp_GRCh38, "Low coverage regions GRCh38", col="black")

dev.off()


