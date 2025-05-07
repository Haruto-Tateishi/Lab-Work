# Synonymous Codon

1. Overview
	- Within the variants on vcf that are annotated as synonymous mutation, generate SFS for each synonymous codon pair (e.g. tgC->tgT) and assess the difference in the strength of natural selection between synonymous codon pairs.
	
2. Dataset
	- 1kG chr1-22(GRCh37.75)
	
3. Annotation Tools 
	- VEP
	
4. Pipeline
	-  Annotate VCF files
		- By using the annotation tools, annotate VCF files about the consequences (e.g. synonymous variant, intron variant)
		- Example usage:
			- ./vep -i path/to/input.vcf --cache --assembly GRCh37 --vcf --pick --pickorder canonical,tsl -o path/to/output.vcf
	- Filter out non-ancestral variants
		- By using chimp FASTA file and a script, vcf_filter_chimp.py, keep the variants whose reference base matches the reference variant of chimpanzee.
	- Build a super dictionary
		- Run vcf_pickle.py
		- Super Dictionary: 1st key; Synonymous Codon Pair (e.g. cgC->cgT), 2nd key; AC_AN, value; counts
		- VEP annotation for synonymous variants contain synonymous codon pair information. However, you can run vcf_ann_codon.py, which annotates flanking bases of the variants based on FASTA file. 
	- Generate SFS
		- Run vcf_pickle_SFS.py
		- It will generate SFS for each consequence. 
	- Visualize relative cumulative SFS (Optional)
		- Run a supplemental script, SFS_rel_cum.py
		- You will obtain relative cumulative graph of SFS.
	- Integration of relative cumulative SFS (Optional)
		- Run a supplemental script, SFS_integration.py
		- You will obtain area under curves of relative cumulative SFS.