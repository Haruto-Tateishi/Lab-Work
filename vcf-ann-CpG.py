# this script adds infromation of CpG status on vcf file.
# you can change the order of annotation by changing the index of insert

import pysam

fn1 = "1kG-chr22-hg19-snpEff-prcsd-codon-reg-exon.vcf"
f1 = pysam.VariantFile(fn1)

fn2 = "Chr22-CpG-list.tsv"
f2 = open(fn2, "r")

if 'CpG' not in f1.header.info:
        f1.header.info.add('CpG', 1, 'String', 'CpG Island Status')

fn3 = "1kG-chr22-hg19-snpEff-prcsd-codon-reg-exon-CpG.vcf"
f3 = pysam.VariantFile(fn3, "w", header=f1.header)

CpG_lines = f2.readlines()
i = 0

for record in f1:
# for g in range(1000):
  vcf_pos = record.pos
  print(vcf_pos)
  while True:
    try:
      CpG_row = CpG_lines[i].split(sep="\t")
      # print(reg_row)
      if CpG_row[0] == "22":
        CpG_pos = int(CpG_row[1])
        if CpG_pos < vcf_pos:
          i += 1
          continue
        if CpG_pos > vcf_pos:
          # no CpG SNPs
          record.info['CpG'] = 'F'
          f3.write(record)
          break
        if CpG_pos == vcf_pos:
          # CpG SNPs
          record.info['CpG'] = 'T'
          f3.write(record)
          i += 1
          break
    except IndexError:
      # rest of no CpG SNPs
      record.info['CpG'] = 'F'
      f3.write(record)
      break