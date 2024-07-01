# this script annotates file with regulatory region.

import pysam

fn1 = "1kG-chr22-hg19-snpEff-prcsd-codon.vcf"
f1 = pysam.VariantFile(fn1)

fn2 = "Chr22-reg-list.tsv"
f2 = open(fn2, "r")

if 'Reg' not in f1.header.info:
        f1.header.info.add('Reg', 1, 'String', 'True/False for Regulatory Region')
if 'Src' not in f1.header.info:
        f1.header.info.add('Src', 1, 'String', 'Sources of Regulatory Region')

fn3 = "1kG-chr22-hg19-snpEff-prcsd-codon-reg.vcf"
f3 = pysam.VariantFile(fn3, 'w', header=f1.header)

reg_lines = f2.readlines()
i = 0

for record in f1:
# for g in range(1000):
  vcf_pos = record.pos
  print(vcf_pos)
  while True:
    try:
      reg_row = reg_lines[i].split(sep="\t")
      # print(reg_row)
      if reg_row[0] == "22":
        reg_pos = int(reg_row[1])
        if reg_pos < vcf_pos:
          i += 1
          continue
        if reg_pos > vcf_pos:
          # no regulatory-region SNPs
          record.info['Reg'] = "F"
          record.info['Src'] = 'N/A'
          f3.write(record)
          break
        if reg_pos == vcf_pos:
          # regulatory-region SNPs
          record.info['Reg'] = 'T'
          record.info['Src'] = reg_row[2]
          i += 1
          f3.write(record)
          break
    except IndexError:
      # rest of no regulatory-region SNPs
      record.info['Reg'] = 'F'
      record.info['Src'] = 'N/A'
      f3.write(record)
      break

f1.close()
f2.close()  
f3.close()