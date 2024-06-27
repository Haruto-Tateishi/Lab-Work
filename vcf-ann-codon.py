# this script add info about codon sequence

import pysam
from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from Bio.Data import CodonTable
import pandas as pd

gtf_file = 'hg19.gtf'
gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, names=[
    'seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

def get_strand(chrom, pos):
    exon = gtf[(gtf['seqname'] == chrom) & (gtf['feature'] == 'exon') & (gtf['start'] <= pos) & (gtf['end'] >= pos)]
    # print(exon)
    if not exon.empty:
        # print(exon.iloc[0]['strand'], exon.iloc[0]['start'], exon.iloc[0]['end'])
        return exon.iloc[0]['strand'], exon.iloc[0]['start'], exon.iloc[0]['end']
    else:
        # print(chrom, pos)
        return None, None, None

reference = pysam.FastaFile('hg19.fa')

vcf_in = pysam.VariantFile('1kG-chr22-hg19-snpEff-prcsd.vcf')
contigs = reference.references

for contig in contigs:
    vcf_in.header.contigs.add(contig, length=reference.get_reference_length(contig))

for field in ['RefSeq', 'MutSeq', 'RefAmi', 'MutAmi']:
    if field not in vcf_in.header.info:
        vcf_in.header.info.add(field, 1, 'String', f'{field} from reference and alternate codons')

vcf_out = pysam.VariantFile('1kG-chr22-hg19-snpEff-prcsd-codon.vcf', 'w', header=vcf_in.header)


def get_codon_annotation(variant, ref_seq, strand, cds_start, cds_end):
    contig = 'chr' + variant.chrom
    
    if strand == '+':
        offset = (record.pos - cds_start) % 3
        region_seq = ref_seq.fetch(contig, record.pos - offset - 1, record.pos + 2 - offset)
    if strand == '-':
        offset = (cds_end - record.pos) % 3
        region_seq = ref_seq.fetch(contig, record.pos - offset - 1, record.pos + 2 - offset)
        region_list = reversed(list(region_seq)) 
        region_seq = "".join(region_list)
    if strand == None:
        region_seq = ref_seq.fetch(contig, record.pos - 2, record.pos + 1)
        offset = 1
    
    ref_base = record.ref.upper()
    alt_bases = variant.alts[0].split(',')
    
    if region_seq[offset].upper() != ref_base:
        print(f"Reference base mismatch at {contig}:{record.pos}. Expected {ref_base}, got {region_seq[1]}")
        return None, None, None, None
    else:
        if region_seq[offset].isupper():
            ref_up = True
        if region_seq[offset].islower():
            ref_up = False

    for alt_base in alt_bases:
        if len(alt_base) == len(ref_base):
            mut_codon = list(region_seq)
            if ref_up:
                mut_codon[offset] = alt_base
                mut_codon = ''.join(mut_codon)
                break
            else:
                mut_codon[offset] = alt_base.lower()
                mut_codon = ''.join(mut_codon)
                break

    ref_codon = region_seq

    ref_aa = str(Seq(ref_codon).translate())
    mut_aa = str(Seq(mut_codon).translate())
    
    return ref_codon, mut_codon, ref_aa, mut_aa


for record in vcf_in:
    if len(record.ref) != 1 or any(len(alt) != 1 for alt in record.alts):
        vcf_out.write(record)
        continue
    contig = 'chr' + record.chrom
    strand, cds_start, cds_end = get_strand(record.chrom, record.pos)
    # print(strand, cds_start, cds_end)
    ref_seq, mut_seq, ref_aa, mut_aa = get_codon_annotation(record, reference, strand, cds_start, cds_end)
    
    if not ref_seq or not mut_seq:
        continue
    
    record.info['RefSeq'] = ref_seq
    record.info['MutSeq'] = mut_seq
    record.info['RefAmi'] = ref_aa
    record.info['MutAmi'] = mut_aa

    vcf_out.write(record)

vcf_out.close()