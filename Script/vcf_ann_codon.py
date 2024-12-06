# this script add info about codon sequence

import pysam
from Bio.Seq import Seq
import gzip
# import pandas as pd

# gtf_file = 'hg19.gtf'
# gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, names=[
#     'seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

# def get_strand(chrom, pos):
#     exon = gtf[(gtf['seqname'] == chrom) & (gtf['feature'] == 'exon') & (gtf['start'] <= pos) & (gtf['end'] >= pos)]
#     # print(exon)
#     if not exon.empty:
#         # print(exon.iloc[0]['strand'], exon.iloc[0]['start'], exon.iloc[0]['end'])
#         return exon.iloc[0]['strand'], exon.iloc[0]['start'], exon.iloc[0]['end']
#     else:
#         # print(chrom, pos)
#         return None, None, None


def get_codon_annotation(variant, ref_seq):
    contig = 'chr' + variant.chrom
    region_seq = ref_seq.fetch(contig, variant.pos - 2, variant.pos + 1)
 
    # print(record.pos, region_seq)

    ref_base = variant.ref.upper()
    # alt_bases = variant.alts[0].split(',')
    
    if region_seq[1].upper() != ref_base:
        print(f"Reference base mismatch at {contig}:{variant.pos}. Expected {ref_base}, got {region_seq[1]}")
        return None
    # else:
    #     if region_seq[1].isupper():
    #         ref_up = True
    #     if region_seq[1].islower():
    #         ref_up = False

    # for alt_base in alt_bases:
    #     if len(alt_base) == len(ref_base):
    #         mut_codon = list(region_seq)
    #         if ref_up:
    #             mut_codon[1] = alt_base
    #             mut_codon = ''.join(mut_codon)
    #             break
    #         else:
    #             mut_codon[1] = alt_base.lower()
    #             mut_codon = ''.join(mut_codon)
    #             break

    ref_codon = region_seq.upper()
    # print(ref_codon)

    # ref_aa = str(Seq(ref_codon).translate())
    # mut_aa = str(Seq(mut_codon).translate())

    return ref_codon


def codon_annotation(reference):
    contigs = reference.references
    for chr in range(1,22):
        value_skipped = 0
        os_skipped = 0
        format_skipped = 0
        non_seq_skipped = 0
        gunzipped = gzip.open(f'vcf/1kG_chr{chr}_vep_pickorder_pantro6filtered_regfiltered_modified.vcf.gz', 'rb')
        vcf_in = pysam.VariantFile(gunzipped)
        # for contig in contigs:
        #     vcf_in.header.contigs.add(contig, length=reference.get_reference_length(contig))
        # for field in ['RefSeq']:
        #     if field not in vcf_in.header.info:
        #         vcf_in.header.info.add(field, 1, 'String', f'{field} from reference and alternate codons')
        vcf_out = pysam.VariantFile(f'vcf/1kG_chr{chr}_vep_pickorder_pantro6filtered_regfiltered_modified_codon.vcf.gz', 'w', header=vcf_in.header)
        while True:
            try:
                record = next(vcf_in)
            except OSError:
                os_skipped +=1
                continue
            except StopIteration:
                break
            # print(record.pos)
            if len(record.ref) != 1 or any(len(alt) != 1 for alt in record.alts):
                # vcf_out.write(record)
                format_skipped+=1
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

codon_annotation(reference)