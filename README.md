# FLASH-Seq

## Information:

This repository contains the scripts related to the analysis of the FLASH-Seq data.

If something is still missing please just contact me !

## Folder Architecture:

```
* pre-processing
* UMI_in_R2
* variant_calling
* R
	- UMI_and_strandInvasion
	- archive_first_draft_to_keep
	- extended_files
	- general
	- timing_graph
	- w18_organoids
* experimental
```

## Folder Description:

* pre-processing
	* Demultiplexing, mapping, counting, 10x genomics, ...
* UMI_in_R2
	* Search of UMI reads in read 2.
* variant_calling
* R/UMI_and_strandInvasion
	* Scripts related to the strand-invasion analysis.
	* **filterInvasionEventsfromBAM.R**: remove suspected strand-invasion events from BAM file.
* R/archive_first_draft_to_keep
	* Legacy code used in the first draft. 
* R/extended_files
	* Scripts related to the development of FLASH-Seq (=extended_file_1).
* R/general
	* FLASH-Seq analysis.
	* Miscellanous functions to process / aggregate files.
* R/timing_graph
	* Arrow-graphs
* R/w18_organoids
	* Post-processing of the week-18 retinal organoids data.
* experimental
	* On-going development, not used in the manuscript.
	* **FS_fastq_to_zUMIs_unmappedBAM.R**: bridge between FLASH-Seq UMI data and zUMIs in order to use the Smart-seq3 transcript reconstruction tool.


