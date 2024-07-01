# this script creates annotation for simple repeats based on the letter size of codon sequence.

import pysam

fn1 = "1kG-chr22-hg19-snpEff-prcsd-codon-reg-exon-CpG.vcf"
f1 = pysam.VariantFile(fn1)

if 'Rep' not in f1.header.info:
        f1.header.info.add('Rep', 1, 'String', 'Simple Repeat Region/Non-Simple Repeat Region')

fn2 = "1kG-chr22-hg19-snpEff-prcsd-codon-reg-exon-CpG-repeat.vcf"
f2 = pysam.VariantFile(fn2, "w", header=f1.header)

for record in f1:
    codon_seq = record.info['RefSeq']
    ref_base = list(codon_seq)[1]
    print(record.pos)
    ref_up = ref_base.isupper()
    if ref_up:
        # non-simple repeat region
        record.info['Rep'] = "F"
        f2.write(record)
        continue
    else:
        # simple repeat region
        record.info['Rep'] = "T"
        f2.write(record)
        continue

          
