# this script add info about codon sequence into the vcf file

import pysam
from Bio.Seq import Seq
import gzip
import pandas as pd
import random

# get the strand and cds start and end positions for a given chromosome and position
def get_codon_annotation(variant, ref_seq):
    contig = 'chr' + variant.chrom
    region_seq = ref_seq.fetch(contig, variant.pos - 2, variant.pos + 1)
 
    ref_base = variant.ref.upper()
       
    if region_seq[1].upper() != ref_base:
        print(f"Reference base mismatch at {contig}:{variant.pos}. Expected {ref_base}, got {region_seq[1]}")
        return None

    ref_codon = region_seq.upper()

    return ref_codon

# create a new VCF file with the codon annotation
def codon_annotation(reference):
    contigs = reference.references
    for chr in range(22,23):
        value_skipped = 0
        os_skipped = 0
        format_skipped = 0
        non_seq_skipped = 0
        zero_one = [0,1]
        gunzipped = gzip.open(f'/mnt/shared_dir/vcf/uk10k/uk10k_pantro6filtered_regfiltered/uk10k_chr{chr}_official_merged_vep_pickorder_fixed_pantro6filtered_regfiltered_modified.vcf.bgz', 'rb')
        vcf_in = pysam.VariantFile(gunzipped)
        for contig in contigs:
            vcf_in.header.contigs.add(contig, length=reference.get_reference_length(contig))
        for field in ['RefSeq']:
            if field not in vcf_in.header.info:
                vcf_in.header.info.add(field, 1, 'String', f'{field} from reference and alternate codons')
        vcf_out = pysam.VariantFile(f'/mnt/shared_dir/vcf/uk10k/uk10k_pantro6filtered_regfiltered/uk10k_pantro6filtered_regfiltered_codon/uk10k_chr{chr}_official_merged_vep_pickorder_fixed_pantro6filtered_regfiltered_modified_codon.vcf', 'w', header=vcf_in.header)
        while True:
            try:
                record = next(vcf_in)
            except OSError:
                os_skipped +=1
                continue
            except StopIteration:
                break
            if len(record.ref) != 1 or any(len(alt) != 1 for alt in record.alts):
                # vcf_out.write(record)
                format_skipped+=1
                print(random.choice(zero_one))
                continue
            contig = 'chr' + record.chrom
            # strand, cds_start, cds_end = get_strand(record.chrom, record.pos)
            # print(strand, cds_start, cds_end)
            ref_seq = get_codon_annotation(record, reference, )
            
            if not ref_seq:
                non_seq_skipped+=1
                continue
            
            record.info['RefSeq'] = ref_seq
            # record.info['MutSeq'] = mut_seq
            # record.info['RefAmi'] = ref_aa
            # record.info['MutAmi'] = mut_aa

            try:
                vcf_out.write(record)
            except ValueError:
                value_skipped +=1
        print(f"{chr} done os_skipped:{os_skipped} value_skipped:{value_skipped} format_skipped:{format_skipped} non_seq_skipped:{non_seq_skipped}")
        vcf_in.close()
        vcf_out.close()

reference = pysam.FastaFile('Document/Other/hg19.fa')

# main function
codon_annotation(reference)