# This script filters VCF files based on reference sequences from a FASTA file in this case chimp FASTA file.

import pysam
import gzip
from pyfaidx import Fasta
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def filter_vcf(vcf_file, fasta_file, output_file, chrom):
    # Open the FASTA file
    genome = Fasta(fasta_file)

    # Open the input VCF file
    gunzipped = gzip.open(vcf_file, 'r')
    vcf = pysam.VariantFile(gunzipped)

    # Open the output VCF file
    with pysam.BGZFile(output_file, 'w') as out_vcf:
        # Write the header
        out_vcf.write(str(vcf.header).encode())

        # Counter for processed variants
        processed = 0
        matched = 0

        # Iterate through the VCF records
        for record in vcf:
            processed += 1
            if processed % 100000 == 0:
                logging.info(f"Processed {processed} variants") 

            # Get the reference base from the FASTA file
            fasta_base = genome[chrom][record.pos - 1].seq.upper()

            # Check if the VCF reference matches the FASTA reference
            if record.ref.upper() == fasta_base:
                # Write the matching record to the output file
                out_vcf.write(str(record).encode())
                matched +=1

    logging.info(f"Finished processing. Total variants processed: {processed}. Matched: {matched}")

# File paths
fasta_file = "pantro6_hg19.fa"



for i in range(1, 23):
    vcf_file = "/home/haruto/ensembl-vep/output.vcf/1kG-chr{}-vep-every.vcf.gz".format(i)
    output_file = "1kG-chr{}-pantro6-filtered.vcf.bgz".format(i)
    chrom = "chr{}".format(i)

    # Run the filtering
    # logging.info("Starting VCF filtering")
    filter_vcf(vcf_file, fasta_file, output_file, chrom)
    # logging.info("VCF filtering complete")

    # Index the output file
    logging.info("Indexing the output VCF file")
    pysam.tabix_index(output_file, preset="vcf", force=True)
    logging.info("Indexing complete")
