# Mutation Rate

1. Overview 
	- Generate SFS by pairing up synonymous variants and intergenic variants, this time, based on mutation rate.
	
2. Dataset
	- 1kG chr1-22(GRCh37.75)
	- roulette (http://genetics.bwh.harvard.edu/downloads/Vova/Roulette/hg19/)
	
3. Annotation Tools 
	- VEP

4. Mask File
	- Pilot
		- 20141020.pilot_mask.whole_genome.bb
	- Strict
		- 20141020.strict_mask.whole_genome.bb
	
5. Pipeline
	- Annotate VCF files
		- By using the annotation tools, annotate VCF files about the consequences (e.g. synonymous variant, intron variant)
		- Example usage:
			- ./vep -i path/to/input.vcf --cache --assembly GRCh37 --vcf --pick --pickorder canonical,tsl -o path/to/output.vcf
	- Filter out non-ancestral variants
		- By using chimp FASTA file and a script, vcf_filter_chimp.py, keep the variants whose reference base matches the reference variant of chimpanzee.
	- Mutation rate
		- For synonymous variant in the filtered vcf, generate a super dictionary which contains mutation rate from roulette vcf file.
	- Bin number for mutation rate
		- Generate bin eedges so that each bin has similar number of variants
	- Build a new super dictionary
		- Pair up selected and neutral variants based on position, bin number, Synonymous Codon Pair/GC Bias.
			- GC Bias has three categories, A:A/T -> C/G, B:C/G -> A/T, n:C/G -> C/G or A/T -> A/T
		- Two ways to pair them up
			- One by One: Pick up the closest neutral for each selected variant one by one.
			- Closest: For the entire list of positions, pair up the variants first from smallest distances to the furthest distances.
				- Closest option is generally better. Thus, it was used for analysis of mask-file-filtered data for both GC Bias and Codon Pair. 
	- Generate SFS
		- Run vcf_pickle_SFS.py
		- It will generate SFS for each GC Bias or Synonymous Codon Pair. 
	- Visualize relative cumulative SFS (Optional)
		- Run a supplemental script, SFS_rel_cum.py
		- You will obtain relative cumulative graph of SFS.
	- Integration of relative cumulative SFS (Optional)
		- Run a supplemental script, SFS_integration.py
		- You will obtain area under curves of relative cumulative SFS.
