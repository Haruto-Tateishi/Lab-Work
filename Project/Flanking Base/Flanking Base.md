# Note of Project

## Project with Jody
1. Overview of the Project with Jody
	- Generate SFS by pairing up synonymous variants and intergenic variants (selected and neutral) based on flanking bases. 
	- Test SF Ratios of each SFS 
		
2. Dataset
	- uk10k chr1-22(GRCh37.75)
	- 1kG chr1-22(GRCh37.75)
	
3. Annotation Tools 
	- VEP
	
4. BED file (Downloaded from UCSC genome browser)
	- ENCODE TFBS
	- FANTOM5
	- OrgAnno
		
5. SF Ratios
	- https://github.com/jodyhey/SF_Ratios

6. Pipeline
	- Annotate VCF files
		- By using the annotation tools, annotate VCF files about the consequences (e.g. synonymous variant, intron variant)
		- Example usage:
			- ./vep -i path/to/input.vcf --cache --assembly GRCh37 --vcf --pick --pickorder canonical,tsl -o path/to/output.vcf
	- Filter out non-ancestral variants
		- By using chimp FASTA file and a script, vcf_filter_chimp.py, keep the variants whose reference base matches the reference variant of chimpanzee.
	- Filter out regulatory region
		- Use bedtool and combine 3 bed files and run vcf_bed.py.
		- Filter out the variants whose positions fall into the intervals of combined BED file.
	- Annotate flanking bases
		- Run vcf_ann_codon.py 
		- It will provide the nucleotide adjacent to the reference base. 
	- Build a super dictionary
		- Run vcf_pickle.py
		- Super dictionary: 1st key; Flanking base, 2nd key; Consequence, 3rd key; position, value; Allele counts.
		- In this case you will use two consequences, synonymous and intergenic variants.
	- Generate SFS
		- Run vcf_pickle_SFS.py
		- It will generate SFS for each consequence by pairing up closest selected variants and neutral variants. 
	- SF Ratios
		- Run SF_Ratios_run.py
		- You will get output data for SF Ratios
		- You can compare those files by running SF_Ratios_analyze.py
	- Visualize relative cumulative SFS (Optional)
		- Run a supplemental script, SFS_rel_cum.py
		- You will obtain relative cumulative graph of SFS.
	- Integration of relative cumulative SFS (Optional)
		- Run a supplemental script, SFS_integration.py
		- You will obtain area under curves of relative cumulative SFS.
		