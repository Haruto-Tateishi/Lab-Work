#  Projects

## Generate SFS and Elucidate the Natural Selection on Synonymous Mutation.

1. Consequence
	- By using annotation tools for vcf files, such as VEP and Sift4g, generate SFS for each mutational consequence including synonymous variant and explore the natural selection. 
	
2. Synonymous Codon
	- Within the variants on vcf that are annotated as synonymous mutation, generate SFS for each synonymous codon pair (e.g. tgC->tgT) and assess the difference in the strength of natural selection between synonymous codon pairs.

3. Flanking Base
	- Generate SFS by pairing up synonymous variants and intergenic variants (selected and neutral) based on flanking bases. 
	- Test SF Ratios of each SFS 

4. Mutation Rate
	- Generate SFS by pairing up synonymous variants and intergenic variants, this time, based on mutation rate.
	
## General Materials

### VCF files

1. 1kG 
	-hg19
	-https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/
	
2. uk10k
	- hg19
	- https://www.uk10k.org/

3. GnomAD
	- hg38
	- https://gnomad.broadinstitute.org/

### Annotation Tool

1. Ensembl Variant Effect Predictor (VEP)
	- https://useast.ensembl.org/info/docs/tools/vep/index.html
	