# this script annotates file with regulatory region.

fn1 = "1kG.chr22.SNPEff.prcsd.vcf"
f1 = open(fn1, "r")

fn2 = "Chr22.reg.list.tsv"
f2 = open(fn2, "r")

fn3 = "1kG.chr22.SNPEff.prcsd.reg.vcf"
f3 = open(fn3, "w")

reg_lines = f2.readlines()
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
        reg_row = reg_lines[i].split(sep="\t")
        # print(reg_row)
        if reg_row[0] == "22":
          reg_pos = int(reg_row[1])
          if reg_pos < vcf_pos:
            i += 1
            continue
          if reg_pos > vcf_pos:
            # no regulatory-region SNPs
            F_reg_info = vcf_row[7].split(sep="|")
            F_reg_info.insert(5, "REG=F")
            F_reg_info.insert(6, "N/A")
            vcf_row[7] = "|".join(F_reg_info)
            vcf_F_reg = "\t".join(vcf_row)
            f3.write(vcf_F_reg)
            break
          if reg_pos == vcf_pos:
            # regulatory-region SNPs
            T_reg_info = vcf_row[7].split(sep="|")
            T_reg_info.insert(5, "REG=T")
            T_reg_info.insert(6, reg_row[2])
            vcf_row[7] = "|".join(T_reg_info)
            vcf_T_reg = "\t".join(vcf_row)
            f3.write(vcf_T_reg)
            i += 1
            break
      except IndexError:
        # rest of no regulatory-region SNPs
        F_reg_info = vcf_row[7].split(sep="|")
        F_reg_info.insert(5, "REG=F")
        F_reg_info.insert(6, "N/A")
        vcf_row[7] = "|".join(F_reg_info)
        vcf_F_reg = "\t".join(vcf_row)
        f3.write(vcf_F_reg)
        break