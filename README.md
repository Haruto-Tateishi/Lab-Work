# Note of Project

## Project with Jody
1. Overview of the Project with Jody
	- Get SFSs of gnomad vcf files without downsampling them, SFSs of synonymous codon pairs from 1kG vcf files, and simulated SFSs for multiple values of g.
	
2. Dataset
	- gnomad chr1-22(GRCh38)
	- 1kG chr1-22(GRCh37.75)
	
3. Annotation Tools 
	- VEP(Variant Effective Predictor)
	
4. Explanation
	- simulated SFSs
		- Simulated for different g values, -100, -10, -1, 0, 1, 10, and 100.
		- Tried it for different population size, 1000, 10000, and 100000.
		- Put theta value as it's 1/50 of population size and maxi as it's 9/10 of population size.
		- Use None for misspec and True for returnexpected.
		- All SFS graphs(Cum_sim_~) is available in the file, Charts, on github archive.
	
	- Gnomad SFSs
		- Downsampled for the nc values, 5008, 10000, and 1000000.
		- All SFS graphs(Cum_gnomad_~) is available in the file, Charts, on the github archive.
		- All SFS text files(gnomad-vep-downsmaple-~) is available on the archive.
		- All tables of integration(gnomad-SFS-downsample-integrate-~) is available on the archive. The tables of area under curves are decreasing order, as scrolling down. 
		
	- 1kG SFSs
		- Made an SFS for all the data, comparing each consequence annotated by vep.
		- Filtered the SNPs whose reference bases are not same as the ancestral bases(chimpanzee:pantro6).
		- Filtered the SNPs which do not have 'synonymous_variant' annotation.
		- Made two pickle files, "all" is using all filtered data, "single" is using only lines having one annotation, which is 'synonymous_variant'
		- Made SFSs from the pickle files and from downsampled data with nc = 30, 50, 100, 200, 500.
		- All SFS graphs(Cum_1kG_~) is available in the file, Charts, on the github archive.
		- All SFS text files(1kG-vep-syn-all/single-codon-downsmaple-~) is available on the archive.
		- All tables of integration(1kG-SFS-downsample-integrate-syn~) is available on the archive. The tables of area under curves are decreasing order, as scrolling down. 

5. Scripts
	- 1kG-pickle.py
		- Create super dictionary, where the first keys are codon pairs(changeable), and the second keys are counts of AN and AC and the values are the counts of appearance, and save it as pickle file.
	- SFS-integration.py
		- Make table file of area under curves for each first key of super dictionary from pickle file or each SFS.
		- It's decreasing order as it goes down to the bottom.
		- Able to integrate simulated SFS
 	- SFS-rel-cum.py
 		- Make visual graph of relative cumulative SFS from pickle file or SFS text file. 
 	- SFS-sim.py
 		- Make visual graph of relative cumulative SFS from simulated SFS.
 	- table-combine.py
 		- Create another table by combining the tables of area under curve for further comparison.
 	- vcf-SFS-conseq-type.py
 		- Make various type of SFS such as bar charts and scatter charts
 	- vcf-ann-codon.py
 		- Add equence information around the reference bases in the info field of vcf file.
 	- vcf-dict-pickle.py
 		- Used for creating super-dictionary pickle file based with consequences as the first keys.
 	- vcf-downsample-pickle.py
 		- From super-dictionary pickle file, create SFS text files by downsampling the data to specific nc value.
 	- vcf-filter-chimp.py
 		- Filter out lines whose ref bases do no match to ones of pantro6 ancestral genome sequence.
 	- vcf-filter-syn.py
 		- Filter out lines which do not have "synonymous_variant" vep annotation.
 	- vcf-snpEff-VEP-types.py
 		- Extract all types of consequences annotated by either VEP or snpEff and output as text file.
 
 6. Documents
 	- Pickle
 		- Contain pickle files.
 		- Could not upload pickle file from gnomad vcf files due to the huge size of the file.
 	- Table
 		- Contain .tsv files such as the tables of curves integration.
 	- SFS
 		- Contain SFS text file.
 	- Other
 		- Contain additional files such as synonymous_codon_pairs.txt
 		
 7. Pipeline
 	- Getting SFS for consequences
 		1. Annotate vcf files with VEP(variant effective predictor) like the command below
 			- ./vep -i path/to/input_vcf --cache --assembly GRCh37 --vcf --everything -o path/to/output_vcf
 		2. Run the script, 1kG-pickle.py or vcf-dict-pickle.py, to get super-dictionary pickle file.
 		3. Run the script, SFS-rel-cum.py, to get visual relative cumulative SFS.
 		4. Run the script, SFS-integration.py, to get table of area under curve for the comparison of degrees of selection.
 	- Getting SFS for synonymous codon pairs
 		1. Annotate vcf files with VEP(variant effective predictor) like the command below
 			- ./vep -i path/to/input_vcf --cache --assembly GRCh37 --vcf --everything -o path/to/output_vcf
 		2. Run the script, vcf-filter-chimp.py, to filter out lines whose ref bases do no match to ones of pantro6 ancestral genome sequence
 		3. Run the script, vcf-filter-syn.py, to filter out lines which do not have "synonymous_variant" vep annotation
 		4. Run the script, 1kG-pickle, to get super-dictionary pickle file.
 		5. Run the script, vcf-downsample-pickle.py, with optional nc values to downsample the pickle file.
 		6. Run the script, SFS-rel-cum, to get visual relative cumulative SFS for each SFS text files and the original pickle file.
 		7. Run the script, SFS-integration.py, to get tables of area under curve for the comparison of degrees of selection.
 		8. Run the script, table-combine.py, to create another table by combining the tables of area under curve for further comparison.
 			
 			
## Project with Vitor		
1. Overview of the Project with Vitor
	- Infer human (non)synonymous DFE(s) using the PRFratio method. For this, we need to define a set of SNPs potentially neutral. In Drosophila melanogaster, we saw that SNPs from short-introns (< 86 bp) were a good candidate. For Humans, we are going to start testing SNP from non-regulatory, intergenic regions and use the folded SFS derived from these SNPs as the denominator for the ratio between two SFSs, where the numerator would be either a synonymous or a nonsynonymous SFS. The goal is to obtain a better version of DFEs that accounts for demography (a characteristic of the PRFratio method) and that uses a more appropriate set of neutral SNPs. It opens for the possibility of measuring selection on different genomic features like: regulatory regions, transcription factor binding sites, UTRs etc, and for measuring selection on synonymous codon pairs (one of the goals of the multiclass synonymous sites project). 
	- All scripts, toy vcf files, and other files is available in the file, Vitor, on the archive.
	
2. BED files I used
	- geneHancer(reg)
	- ORegAnno(reg)
	- ENCODE(reg)
	- NCBI RefSeq(exon/intron)
	- CpG island(CpG)
	- FASTA file for hg19

3. Human Genome ID
	- GRCh37.75/hg19

4. Annotation tools for SNPs (synonymous, nonsynonymous, intergenic, etc)
	- SNPEff

5. Dataset I used
	- vcf file from 1000 genome project, phase 3, hg19
	- 5 BED files and other files from UCSC genome browser

6. Sample information of 1kG vcf
	- super population
		- population
	- African (AFR)
		- YRI	Yoruba	Yoruba in Ibadan, Nigeria
		- LWK	Luhya	Luhya in Webuye, Kenya
		- GWD	Gambian	Gambian in Western Division, The Gambia
		- MSL	Mende	Mende in Sierra Leone
		- ESN	Esan	Esan in Nigeria
	- Admixed American (AMR)
		- ASW	African-American SW	African Ancestry in Southwest US
		- ACB	African-Caribbean	African Caribbean in Barbados
		- MXL	Mexican-American	Mexican Ancestry in Los Angeles, California
		- PUR	Puerto Rican	Puerto Rican in Puerto Rico
		- CLM	Colombian	Colombian in Medellin, Colombia
		- PEL	Peruvian	Peruvian in Lima, Peru
	- East Asian (EAS)
		- CHB	Han Chinese	Han Chinese in Beijing, China
		- JPT	Japanese	Japanese in Tokyo, Japan
		- CHS	Southern Han Chinese	Han Chinese South
		- CDX	Dai Chinese	Chinese Dai in Xishuangbanna, China
		- KHV	Kinh Vietnamese	Kinh in Ho Chi Minh City, Vietnam
		- CHD	Denver Chinese	Chinese in Denver, Colorado (pilot 3 only)
	- Europeans (EUR)
		- CEU	CEPH	Utah residents (CEPH) with Northern and Western European ancestry
		- TSI	Tuscan	Toscani in Italia
		- GBR	British	British in England and Scotland
		- FIN	Finnish	Finnish in Finland
		- IBS	Spanish	Iberian populations in Spain
	- South Asians (SAS)
		- GIH	Gujarati	Gujarati Indian in Houston, TX
		- PJL	Punjabi	Punjabi in Lahore, Pakistan
		- BEB	Bengali	Bengali in Bangladesh
		- STU	Sri Lankan	Sri Lankan Tamil in the UK
		- ITU	Indian	Indian Telugu in the UK

7. Meeting logs
	- 05/22/24
		- Regulatory region VCF annotation that is done; 
		- Discussed creating a table or getting the SFSs from the VCF; 
		- Discussed the sequence category and issues with closest SNPs; 
		- Discussed using samtools or pysam for creating sequences categories; 
		- Discussed what information the table should have; 
		- Which population should we analyze (and superpopulations)? 
		- Haruto will check the sample size of each of the above populations; 
		- Should we apply our method to the population or subpopulation level? 

8. To-Do list
	- [x] Create table for SFS after all annotations
	- [ ] Create several SFS by changing the filter options