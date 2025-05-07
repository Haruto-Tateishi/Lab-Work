# Consequence

1. Overview
	-  By using annotation tools for vcf files, such as VEP and SIFT4G, generate SFS for each mutational consequence including synonymous variant and explore the natural selection. 
		
2. Dataset
	- gnomad chr1-22(GRCh38)
	- 1kG chr1-22(GRCh37.75)
	
3. Annotation Tools 
	- VEP
	- SIFT4G
		- https://sift.bii.a-star.edu.sg/sift4g/AboutSIFT4G.html
	
4. Pipline
	- Annotate VCF files
		- By using the annotation tools, annotate VCF files about the consequences (e.g. synonymous variant, intron variant)
		- Example usage:
			- ./vep -i path/to/input.vcf --cache --assembly GRCh37 --vcf --pick --pickorder canonical,tsl -o path/to/output.vcf
	- Filter out non-ancestral variants
		- By using chimp FASTA file and a script, vcf_filter_chimp.py, keep the variants whose reference base matches the reference variant of chimpanzee.
	- Build a super dictionary
		- Run vcf_pickle.py
		- Super dictionary: 1st key; Consequence, 2nd key; AC_AN, value; counts.
		- By using a function, conseq_pickle_specific_conseq, you can select specific consequences when you build a super dictionary.
	- Generate SFS
		- Run vcf_pickle_SFS.py
		- It will generate SFS for each consequence. 
	- Visualize relative cumulative SFS (Optional)
		- Run a supplemental script, SFS_rel_cum.py
		- You will obtain relative cumulative graph of SFS.
	- Integration of relative cumulative SFS (Optional)
		- Run a supplemental script, SFS_integration.py
		- You will obtain area under curves of relative cumulative SFS.

