
## Tools

* 1_get_GC_content_from_GTF.R
	* ≈ https://github.com/dpryan79/Answers/blob/master/SEQanswers_42420/GTF2LengthGC.R
* 2_get_exon_intron_from_GTF.R
	* ≈ https://www.biostars.org/p/165226/
	* Extract exon / intron position per gene and reduce/merge overlapping regions to create a new GTF.
	* Used in the mapping of FS-UMI/SS3 to explore the orientation of the reads compared to the reference (exon or intron). 

## Functions

* 4_aggregate_files_functions.R
	* Functions to aggregate various input (featurecounts, STAR.log, ...)
* 6_singleCell_functions.R
	* Functions to process scRNA data. 
* shinyUMAP.R
	* Given a seurat object processed with a standard pipeline and a facultative differential expression (seurat) table, create an interractive shiny window where the user can visualize the seurat cluster, gene expression and differential expression results together. Kinda like UCSC cell browser but with the DEG results and from the R session. Tested with up-to 100K cells. 

## Prepare files

* 5_aggregate_preprocess_results.R
	* Uses 4_aggregate_files_functions.R to aggregate all the files.

## Analysis

* 3_hPBMCs_analysis.Rmd
* 8_HEK_analysis.Rmd
* 7_chr12_coverage_karyoplotr.R
