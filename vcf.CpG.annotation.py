# this script adds infromation of CpG status on vcf file which has already been annotated for SNPs type and reglatory region.
# you can change the order of annotation by changing the index of insert

fn1 = "1kG.chr22.SNPEff.prcsd.reg.vcf"
f1 = open(fn1, "r")

fn2 = "Chr22.CpG.list.tsv"
f2 = open(fn2, "r")

fn3 = "1kG.chr22.SNPEff.prcsd.reg.CpG.vcf"
f3 = open(fn3, "w")

CpG_lines = f2.readlines()
i = 0

while True:
# for g in range(1000):
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
        CpG_row = CpG_lines[i].split(sep="\t")
        # print(reg_row)
        if CpG_row[0] == "22":
          CpG_pos = int(CpG_row[1])
          if CpG_pos < vcf_pos:
            i += 1
            continue
          if CpG_pos > vcf_pos:
            # no CpG SNPs
            F_CpG_info = vcf_row[7].split(sep="|")
            F_CpG_info.insert(7, "CpG=F")
            vcf_row[7] = "|".join(F_CpG_info)
            vcf_F_CpG = "\t".join(vcf_row)
            f3.write(vcf_F_CpG)
            break
          if CpG_pos == vcf_pos:
            # CpG SNPs
            T_CpG_info = vcf_row[7].split(sep="|")
            T_CpG_info.insert(7, "CpG=T")
            vcf_row[7] = "|".join(T_CpG_info)
            vcf_T_CpG = "\t".join(vcf_row)
            f3.write(vcf_T_CpG)
            i += 1
            break
      except IndexError:
        # rest of no CpG SNPs
        F_CpG_info = vcf_row[7].split(sep="|")
        F_CpG_info.insert(7, "CpG=F")
        vcf_row[7] = "|".join(F_CpG_info)
        vcf_F_CpG = "\t".join(vcf_row)
        f3.write(vcf_F_CpG)
        break