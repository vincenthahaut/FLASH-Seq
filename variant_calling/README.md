# FLASH-Seq Variant Calling

* **Mapping / BAM processing / Variant Calling:** scRNA_variant_calling.sh
* **Data Analysis:** scRNA_variant_calling_postprocess.R


Following GATK best practises, bcftools mpileup provided results similar to GATK haplotypecaller and was >10-times faster.

See: https://groups.google.com/g/rna-star/c/Cpsf-_rLK9I to generate the splice-site reference. 

**Notes:** 

* The percentage of variant shared between WES and FS can be further increased by including more WES with a lower coverage (==> often at the start/end of the WES coverage of an exon). 
* The FPR can be decreased (30-50%) by filtering out variants found in <10 cells.

