# this script adds infromation of exon on vcf file.
# you can change the order of annotation by changing the index of insert

import pysam

fn1 = "1kG-chr22-hg19-snpEff-prcsd-codon-reg.vcf"
f1 = pysam.VariantFile(fn1)

fn2 = "Chr22-exon-list.tsv"
f2 = open(fn2, "r")

if 'Exo/Int' not in f1.header.info:
        f1.header.info.add('Exo/Int', 1, 'String', 'Coding Region/Non-Coding Region')
if 'Str' not in f1.header.info:
        f1.header.info.add('Str', 1, 'String', 'Positive Strand/Negative Strand')

fn3 = "1kG-chr22-hg19-snpEff-prcsd-codon-reg-exon.vcf"
f3 = pysam.VariantFile(fn3, "w", header=f1.header)

ex_lines = f2.readlines()
i = 0

for record in f1:
# for g in range(10000):
  vcf_pos = record.pos
  print(vcf_pos)
  while True:
    try:
      ex_row = ex_lines[i].split(sep="\t")
      # print(reg_row)
      if ex_row[0] == "22":
        ex_pos = int(ex_row[1])
        if ex_pos < vcf_pos:
          i += 1
          continue
        if ex_pos > vcf_pos:
          # no exon SNPs, but intron SNPs
          record.info['Exo/Int'] = 'intron'
          record.info['Str'] = "N/A"
          f3.write(record)
          break
        if ex_pos == vcf_pos:
          # exon SNPs
          record.info['Exo/Int'] = 'exon'
          record.info['Str'] = ex_row[2]
          f3.write(record)
          i +=1
          break
    except IndexError:
      # rest of no exon SNPs, but intron SNPs
      record.info['Exo/Int'] = 'intron'
      record.info['Str'] = "N/A"
      f3.write(record)
      break