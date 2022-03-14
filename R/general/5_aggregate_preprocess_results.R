library(tidyverse)
options(scipen=999)
source("/home/vincent.hahaut/Desktop/FLASH-Seq/R/general/4_aggregate_files_functions.R")

# 0. Sequencing Runs
seqPaths <- c(
  "/home/vincent.hahaut/data_storage/200622_NB551561_0024_AH27CVBGXG/STAR/",
  "/home/vincent.hahaut/data_storage/200717_NB551561_0027_AH2HMWBGXG/STAR/",
  "/home/vincent.hahaut/data_storage/200903_NB551561_0029_AH5CWTBGXG/STAR/",
  "/home/vincent.hahaut/data_storage/200904_NB551561_0030_AH5CW5BGXG/STAR/",
  "/home/vincent.hahaut/data_storage/201005_NB551561_0031_AH5CWGBGXG/STAR/",
  "/home/vincent.hahaut/data_storage/201030_NB551561_0034_AHTTVFBGXG/STAR/",
  "/home/vincent.hahaut/data_storage/201113_NB551561_0035_AHTTV7BGXG/STAR/",
  "/home/vincent.hahaut/data_storage/210122_NB551561_0039_AHH5VCBGXH/STAR/",
  "/home/vincent.hahaut/data_storage/210309_NB551561_0043_AHNTN5AFX2/STAR/",
  "/home/vincent.hahaut/data_storage/210322_NB551561_0044_AHNTL2AFX2/STAR/",
  "/home/vincent.hahaut/data_storage/210428_NB551561_0049_AHTYNGAFX2/STAR/",
  "/home/vincent.hahaut/data_storage/210323_NB551561_0045_AHTVN5BGXG/STAR/",
  "/home/vincent.hahaut/data_storage/210617_NB551561_0051_AHTYJYAFX2/STAR/",
  "/home/vincent.hahaut/data_storage/210818_NB551561_0056_AHW77TBGXJ/STAR/",
  "/home/vincent.hahaut/data_storage/210906_NB551561_0058_AH5MLFBGXK/STAR/",
  "/home/vincent.hahaut/data_storage/211014_NB551561_0062_AHF5L2BGXK/STAR/",
  "/home/vincent.hahaut/data_storage/210916_NB551561_0059_AHF5KHBGXK/STAR/",
  "/home/vincent.hahaut/data_storage/211025_NB551561_0066_AHCMHFAFX3/STAR/",
  "/home/vincent.hahaut/data_storage/211019_NB551561_0064_AHCTY5BGXJ/STAR/",
  "/home/vincent.hahaut/data_storage/211110_NB551561_0069_AHMHYLBGXK/STAR_SINGLE/",
  "/home/vincent.hahaut/data_storage/211112_NB551561_0070_AHMHYYBGXK/STAR_SINGLE/",
  "/home/vincent.hahaut/data_storage/220304_NB551561_0086_AHKCVLBGXL/STAR/",
  "/home/vincent.hahaut/data_storage/220303_NB551561_0085_AHK7WTBGXL/STAR/",
  "/home/vincent.hahaut/data_storage/211202_NB551561_0072_AH3YYHBGXK/STAR/"
)

mydate <- Sys.Date()


# 1. STAR Final Logs
paths.STAR.log <- makePath(seqPaths = seqPaths, suffix = "STAR/ID_Log.final.out")
STAR.log <- getSTARLog(paths = paths.STAR.log$paths, sample.ids = paths.STAR.log$sample.ids, threads = 10)
write_rds(STAR.log, file = paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/STAR.logs.", mydate, ".rds.bz"), compress = "bz")

# 2. GC content
paths.GC.log <- makePath(seqPaths = seqPaths, suffix = "ReSQC/ID_GCcontent.txt.GC.xls")
GC.log <- getRESeQC_GC(paths = paths.GC.log$paths, sample.ids = paths.GC.log$sample.ids, threads = 10)
write_rds(GC.log, file = paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/GC.logs.", mydate, ".rds.bz"), compress = "bz")

# 3. ReSQC Coverage
paths.ReSQC.cov <- makePath(seqPaths = seqPaths, suffix = "ReSQC/ID_geneBody.all.geneBodyCoverage.txt")
ReSQC.cov <- getRESeQC_coverage(paths = paths.ReSQC.cov$paths, sample.ids = paths.ReSQC.cov$sample.ids, threads = 10)
write_rds(ReSQC.cov, file = paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/ReSQC.coverage.", mydate, ".rds.bz"), compress = "bz")

# 4. ReSQC Distribution
paths.ReSQC.distr <- makePath(seqPaths = seqPaths, suffix = "ReSQC/ID_readDistribution.txt")
ReSQC.distr <- getRESeQC_readDistribution(paths = paths.ReSQC.distr$paths, sample.ids = paths.ReSQC.distr$sample.ids, threads = 10)
write_rds(ReSQC.distr, file = paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/ReSQC.distribution.", mydate, ".rds.bz"), compress = "bz")

# 7. ReadCounts

# 7.1. ALL - GENES
paths.featurecounts <- makePath(seqPaths = seqPaths, suffix = "FEATURECOUNTS/ID_ReadCount.featureCounts.gencode.GENE.txt")
featurecounts.gene <- getFeatureCount(paths = paths.featurecounts$paths, sample.ids = paths.featurecounts$sample.ids, threads = 10)
write_rds(featurecounts.gene, file = paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/featureCounts.genes.all.", mydate, ".rds.bz"), compress = "bz")

# 7.2. EXONS - INTRON 
paths.featurecounts <- makePath(seqPaths = seqPaths, suffix = "FEATURECOUNTS/ID_ReadCount.featureCounts.gencode.EXONINTRON.DOWNSAMPLING.TIMES.txt", downsampling = 250000)
featurecounts.exonintron <- getFeatureCount(paths = paths.featurecounts$paths, sample.ids = paths.featurecounts$sample.ids, threads = 10)
write_rds(featurecounts.exonintron, file = paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/featureCounts.exonintron.", mydate, ".rds.bz"), compress = "bz")

# 7.2. ALL - EXON 
paths.featurecounts <- makePath(seqPaths = seqPaths, suffix = "FEATURECOUNTS/ID_ReadCount.featureCounts.gencode.txt")
featurecounts.exon <- getFeatureCount(paths = paths.featurecounts$paths, sample.ids = paths.featurecounts$sample.ids, threads = 30)
write_rds(featurecounts.exon, file = paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/featureCounts.exons.all.", mydate, ".rds.bz"), compress = "bz")

# 8. Downsampling

# 8.1. ALL - GENES
mypaths <- makePath(seqPaths = seqPaths, suffix = "FEATURECOUNTS/ID_ReadCount.featureCounts.gencode.GENE.DOWNSAMPLING.TIMES.txt.gz", times = 1, downsampling = c(5000,10000,20000,40000,50000,75000,100000,125000,250000,375000,500000,625000,750000,1000000,1250000,1500000))
ft <- getFeatureCount(paths = mypaths$paths, sample.ids = mypaths$sample.ids, threads = 10)
write_rds(ft, file = paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/downsampling_featureCounts.genes.all.", mydate, ".rds.bz"), compress = "bz")

# 8.2. ALL - EXON
mypaths <- makePath(seqPaths = seqPaths, suffix = "FEATURECOUNTS/ID_ReadCount.featureCounts.gencode.DOWNSAMPLING.TIMES.txt.gz", times = 1, downsampling = c(5000,10000,20000,40000,50000,75000,100000,125000,250000,375000,500000,625000,750000,1000000,1250000,1500000))
ft <- getFeatureCount(paths = mypaths$paths, sample.ids = mypaths$sample.ids, threads = 10)
write_rds(ft, file = paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/downsampling_featureCounts.exon.all.", mydate, ".rds.bz"), compress = "bz")

# 9. Number of expressed genes per cell

# 9.1. ALL - GENES
sce <- read_rds(paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/downsampling_featureCounts.genes.all.", mydate, ".rds.bz"))
med <- sceToExpressedGenes(sce, thresholds = c(0,1,5,10))
write_rds(med, paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/downsampling_medianGenes_featureCounts.genes.all.", mydate, ".rds.bz"))

# 9.2. ALL -EXONS
sce <- read_rds(paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/downsampling_featureCounts.exon.all.", mydate, ".rds.bz"))
med <- sceToExpressedGenes(sce, thresholds = c(0,1,5,10))
write_rds(med, paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/downsampling_medianGenes_featureCounts.exon.all.", mydate, ".rds.bz"))

# 10. Salmon
# I was using a loop but it seem to crash for an unknown reason
options(scipen = 999)

# 10.1. 100K
salmon.path <- makePath(seqPaths = seqPaths, suffix = "salmon_ALL_TRANSCRIPTS_TIMES/quant.sf", times = 1, downsampling = 125000)
salmon <- getSalmon(paths = salmon.path$paths, sample.ids = salmon.path$sample.ids, gtf = "/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.primary_assembly.annotation.gtf", txOUT = TRUE)
write_rds(salmon, file = paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/salmon.all.", 125000, ".", mydate, ".rds.bz"), compress = "bz")

# 10.2. 250K
salmon.path <- makePath(seqPaths = seqPaths, suffix = "salmon_ALL_TRANSCRIPTS_TIMES/quant.sf", times = 1, downsampling = 250000)
salmon <- getSalmon(paths = salmon.path$paths, sample.ids = salmon.path$sample.ids, gtf = "/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.primary_assembly.annotation.gtf", txOUT = TRUE)
write_rds(salmon, file = paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/salmon.all.", 250000, ".", mydate, ".rds.bz"), compress = "bz")

# 10.3. 500K
salmon.path <- makePath(seqPaths = seqPaths, suffix = "salmon_ALL_TRANSCRIPTS_TIMES/quant.sf", times = 1, downsampling = 500000)
salmon <- getSalmon(paths = salmon.path$paths, sample.ids = salmon.path$sample.ids, gtf = "/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.primary_assembly.annotation.gtf", txOUT = TRUE)
write_rds(salmon, file = paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/salmon.all.", 500000, ".", mydate, ".rds.bz"), compress = "bz")

# 11. Number of expressed isoforms per cell

# 11.1. 500K
sce <- read_rds(paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/salmon.all.500000.", mydate, ".rds.bz"))
med <- sceToExpressedGenes(sce, thresholds = c(0,1,5,10,20,50), salmon = TRUE)
write_rds(med, paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/downsampling_medianGenes_salmon.all.500000.", mydate, ".rds.bz"))

# 11.2. 250K
sce <- read_rds(paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/salmon.all.250000.", mydate, ".rds.bz"))
med <- sceToExpressedGenes(sce, thresholds = c(0,1,5,10,20,50), salmon = TRUE)
write_rds(med, paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/downsampling_medianGenes_salmon.all.250000.", mydate, ".rds.bz"))

# 11.3. 125K
sce <- read_rds(paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/salmon.all.125000.", mydate, ".rds.bz"))
med <- sceToExpressedGenes(sce, thresholds = c(0,1,5,10,20,50), salmon = TRUE)
write_rds(med, paste0("/home/vincent.hahaut/data_storage/FS_intermediate_file/downsampling_medianGenes_salmon.all.125000.", mydate, ".rds.bz"))

