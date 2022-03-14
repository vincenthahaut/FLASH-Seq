# FLASH-Seq: Variant Calling

## File Description

* **Mapping / BAM processing / Variant Calling:** scRNA_variant_calling.sh
* **Data Analysis:** scRNA_variant_calling_postprocess.R

## Notes

1. Following GATK best practises, bcftools mpileup provided results similar to GATK haplotypecaller and was >10-times faster.

2. See: https://groups.google.com/g/rna-star/c/Cpsf-_rLK9I to generate the splice-site reference. 

3. The percentage of variant shared between WES and FS can be further increased by including SNPs from WES with a lower quality (==> often at the start/end of the WES coverage of an exon). 

4. The FPR can be decreased by filtering out variants found in <10 cells.

