# this script adds infromation of exon on vcf file which has already been annotated for SNPs type, reglatory region, and CpG status.
# you can change the order of annotation by changing the index of insert

fn1 = "1kG.chr22.SNPEff.prcsd.reg.CpG.vcf"
f1 = open(fn1, "r")

fn2 = "Chr22.exon.list.tsv"
f2 = open(fn2, "r")

fn3 = "1kG.chr22.SNPEff.prcsd.reg.CpG.exon.vcf"
f3 = open(fn3, "w")

ex_lines = f2.readlines()
i = 0

while True:
# for g in range(10000):
  vcf_line = f1.readline()
  if len(vcf_line) ==0:
    break
  if vcf_line[0] == "#":
    f3.write(vcf_line)
    continue
  else:
    vcf_row = vcf_line.split(sep="\t")
    vcf_pos = int(vcf_row[1])
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
            intron_info = vcf_row[7].split(sep="|")
            intron_info.insert(8, "intron")
            intron_info.insert(9, "N/A")
            vcf_row[7] = "|".join(intron_info)
            vcf_intron = "\t".join(vcf_row)
            f3.write(vcf_intron)
            break
          if ex_pos == vcf_pos:
            # exon SNPs
            info = ex_row[2]
            exon_info = vcf_row[7].split(sep="|")
            exon_info.insert(8, "exon")
            exon_info.insert(9, info)
            vcf_row[7] = "|".join(exon_info)
            vcf_exon = "\t".join(vcf_row)
            f3.write(vcf_exon)
            i += 1
            break
      except IndexError:
        # rest of no exon SNPs, but intron SNPs
        intron_info = vcf_row[7].split(sep="|")
        intron_info.insert(8, "intron")
        intron_info.insert(9, "N/A")
        vcf_row[7] = "|".join(intron_info)
        vcf_intron = "\t".join(vcf_row)
        f3.write(vcf_intron)
        break