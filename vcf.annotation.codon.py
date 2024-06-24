# this script add info about codon sequence

import pysam
from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from Bio.Data import CodonTable

reference = pysam.FastaFile('hg38.sequence.fa')

vcf_in = pysam.VariantFile('toy.vcf')
contigs = reference.references

for contig in contigs:
    vcf_in.header.contigs.add(contig, length=reference.get_reference_length(contig))

for field in ['RefSeq', 'MutSeq']:
    if field not in vcf_in.header.info:
        vcf_in.header.info.add(field, 1, 'String', f'{field} from reference and alternate codons')

vcf_out = pysam.VariantFile('annotated_variants.vcf', 'w', header=vcf_in.header)


def get_codon_annotation(variant, ref_seq):
    contig = 'chr' + variant.chrom 
    pos =  variant.pos

    ref_base = ref_seq.fetch(contig, pos, pos + 1)
    print(ref_base)
    print(variant.ref)
    
    if len(variant.alts) == 0 or not all(base in "ACGT" for base in variant.alts[0]):
        return "", "", "", ""
    
    alt_bases = variant.alts[0].split(',')
    # print(alt_bases)
    ref_seq = ref_seq.fetch(contig, pos - 2, pos + 3)
    # print(ref_seq)

    for alt_base in alt_bases:
        if len(alt_base) == len(ref_base):
            mut_seq = list(ref_seq)
            mut_seq[2] = alt_base
            mut_seq = ''.join(mut_seq)
            break
    
    ref_codon_seq = Seq(ref_seq)
    # print(ref_codon_seq)
    mut_codon_seq = Seq(mut_seq)
    # print(mut_codon_seq)
    
    return ref_seq, mut_seq


for record in vcf_in:
    ref_seq, mut_seq = get_codon_annotation(record, reference)
    
    if not ref_seq or not mut_seq:
        continue
    
    # Add annotations to the INFO field of the VCF record
    record.info['RefSeq'] = ref_seq
    record.info['MutSeq'] = mut_seq
    
    vcf_out.write(record)

vcf_out.close()