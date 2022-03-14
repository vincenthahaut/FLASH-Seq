library(karyoploteR)

# 1. Load the Bam files
# These particular BAM were chosen as they depicted well the issue. 
# However, I looked at >100 and almost any other would have been fine as well.
bam.19 <- "/home/vincent.hahaut/data_storage/210323_NB551561_0045_AHTVN5BGXG/STAR/247_HEK_LV1_A6/DOWNSAMPLE/DOWN_500000_1/247_HEK_LV1_A6_downsampled_500000_Aligned.sortedByCoord.filtered.1.bam"
bam.12 <- "/home/vincent.hahaut/data_storage/210323_NB551561_0045_AHTVN5BGXG/STAR/262_HEK_LA_12c_01_C3/DOWNSAMPLE/DOWN_500000_1/262_HEK_LA_12c_01_C3_downsampled_500000_Aligned.sortedByCoord.filtered.1.bam"
bam.10 <- "/home/vincent.hahaut/data_storage/210323_NB551561_0045_AHTVN5BGXG/STAR/261_HEK_LA_10c_00625_E3/DOWNSAMPLE/DOWN_500000_1/261_HEK_LA_10c_00625_E3_downsampled_500000_Aligned.sortedByCoord.filtered.1.bam"
bam.8 <- "/home/vincent.hahaut/data_storage/210323_NB551561_0045_AHTVN5BGXG/STAR/230_HEK_LA_8c_0015625_E4/DOWNSAMPLE/DOWN_500000_1/230_HEK_LA_8c_0015625_E4_downsampled_500000_Aligned.sortedByCoord.filtered.1.bam"
bam.6 <- "/home/vincent.hahaut/data_storage/210323_NB551561_0045_AHTVN5BGXG/STAR/260_HEK_LA_6c_0015625_D3/DOWNSAMPLE/DOWN_500000_1/260_HEK_LA_6c_0015625_D3_downsampled_500000_Aligned.sortedByCoord.filtered.1.bam"
bam.4 <- "/home/vincent.hahaut/data_storage/210323_NB551561_0045_AHTVN5BGXG/STAR/258_HEK_LA_4c_0015625_G3/DOWNSAMPLE/DOWN_500000_1/258_HEK_LA_4c_0015625_G3_downsampled_500000_Aligned.sortedByCoord.filtered.1.bam"

# 2. Create the graph using karyoplotR (chr12)
png("HEK_lowAmplificationCycles_chr12.png", width = 800, height = 600)
pp <- getDefaultPlotParams(plot.type=2)
pp$topmargin <- 150
pp$bottommargin <- 10
pp$data2height <- 1
kp <- plotKaryotype(chromosomes = "chr12", plot.params = pp)
kp <- kpPlotBAMDensity(kp, data=bam.19, col="#B10026", border=NA, r0=0, r1=0.22, window.size = 750000)
kpAxis(kp, ymax=1800, r0=0, r1=0.22, cex=1)
kp <- kpPlotBAMDensity(kp, data=bam.12, col="#E31A1C", border=NA, r0=0.3, r1=0.52, window.size = 750000)
kpAxis(kp, ymax=1800, r0=0.3, r1=0.55, cex=1)
kp <- kpPlotBAMDensity(kp, data=bam.10, col="#FC4E2A", border=NA, r0=0.6, r1=0.82, window.size = 750000)
kpAxis(kp, ymax=1800, r0=0.6, r1=0.82, cex=1)
kp <- kpPlotBAMDensity(kp, data=bam.8, col="#FD8D3C", border=NA, r0=0.9, r1=1.12, window.size = 750000)
kpAxis(kp, ymax=1800, r0=0.9, r1=1.12, cex=1)
kp <- kpPlotBAMDensity(kp, data=bam.6, col="#FEB24C", border=NA, r0=1.2, r1=1.42, window.size = 750000)
kpAxis(kp, ymax=1800, r0=1.2, r1=1.42, cex=1)
kp <- kpPlotBAMDensity(kp, data=bam.4, col= "#FED976", border=NA, r0=1.5, r1=1.72, window.size = 750000)
kpAxis(kp, ymax=1800, r0=1.5, r1=1.72, cex=1)
dev.off()

# 2. Create the graph using karyoplotR (chr12 - centromere)
png("HEK_lowAmplificationCycles_chr12_centromere.png", width = 800, height = 600)
pp <- getDefaultPlotParams(plot.type=1)
pp$topmargin <- 10
pp$bottommargin <- 10
pp$data2height <- 1
pp$data1height <- 750
kp <- plotKaryotype(chromosomes = "chr12", zoom = GRanges("chr12:33200001-37800000"), plot.params = pp)
kp <- kpPlotBAMDensity(kp, data=bam.19, col="#B10026", border=NA, r0=0, r1=0.13, window.size = 50000)
kpAxis(kp, ymax=20, r0=0, r1=0.13, cex=1)
kp <- kpPlotBAMDensity(kp, data=bam.12, col="#E31A1C", border=NA, r0=0.17, r1=0.3, window.size = 50000)
kpAxis(kp, ymax=20, r0=0.17, r1=0.3, cex=1)
kp <- kpPlotBAMDensity(kp, data=bam.10, col="#FC4E2A", border=NA, r0=0.34, r1=0.47, window.size = 50000)
kpAxis(kp, ymax=20, r0=0.34, r1=0.47, cex=1)
kp <- kpPlotBAMDensity(kp, data=bam.8, col="#FD8D3C", border=NA, r0=0.51, r1=0.63, window.size = 50000)
kpAxis(kp, ymax=20, r0=0.51, r1=0.63, cex=1)
kp <- kpPlotBAMDensity(kp, data=bam.6, col="#FEB24C", border=NA, r0=0.67, r1=0.80, window.size = 50000)
kpAxis(kp, ymax=20, r0=0.67, r1=0.80, cex=1)
kp <- kpPlotBAMDensity(kp, data=bam.4, col= "#FED976", border=NA, r0=0.84, r1=0.97, window.size = 50000)
kpAxis(kp, ymax=20, r0=0.84, r1=0.97, cex=1)
dev.off()
