# this script annotates file with regulatory region.

fn1 = "1kG.chr22.SNPEff.prcsd.vcf"
f1 = open(fn1, "r")

fn2 = "1kG.chr22.reg.list.txt"
f2 = open(fn2, "r")

fn3 = "1kG.chr22.SNPEff.prcsd.reg.vcf"
f3 = open(fn3, "w")

reg_lines = f2.readlines()
i = 0

while True:
  vcf_row = f1.readline().split(sep="\t")
  if vcf_row[0] == "":
    break
  if vcf_row[0] == "#CHROM":
    vcf_header = "\t".join(vcf_row)
    f3.write(vcf_header)
    continue
  else:
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
            vcf_row[7] = "|".join(F_reg_info)
            vcf_F_reg = "\t".join(vcf_row)
            f3.write(vcf_F_reg)
            break
          if reg_pos == vcf_pos:
            # regulatory-region SNPs
            T_reg_info = vcf_row[7].split(sep="|")
            T_reg_info.insert(5, "REG=T")
            vcf_row[7] = "|".join(T_reg_info)
            vcf_T_reg = "\t".join(vcf_row)
            f3.write(vcf_T_reg)
            i += 1
            break
      except IndexError:
        # rest of no regulatory-region SNPs
        F_reg_info = vcf_row[7].split(sep="|")
        F_reg_info.insert(5, "REG=F")
        vcf_row[7] = "|".join(F_reg_info)
        vcf_F_reg = "\t".join(vcf_row)
        f3.write(vcf_F_reg)
        break