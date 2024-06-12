# this script makes a list of CpG island SNPs based on BED file
# through running this script all of the numbers of positions will be 1-base

fn1 = "1kG.chr22.SNPEff.prcsd.reg.vcf"
f1 = open(fn1, 'r')

fn2 = "Chr22.CpG.list.tsv"
f2 = open(fn2, 'w')

fn3 = "Chr22.CpG.merged.tsv"
f3 = open(fn3, "r")

BED_lines = f3.readlines()
BED_len = len(BED_lines)
i = BED_len - 1

while True:
# for g in range(1000):
    vcf_data = f1.readline().split(sep="\t")
    if len(vcf_data) == 0:
        break 
    if vcf_data[0] == "22":
        vcf_pos = int(vcf_data[1])
        print(vcf_pos)
        if vcf_pos >= 50800000:
            break
        else:
            while True:
                int_row = BED_lines[i].split(sep="\t")
                # print(int_row)
                if int_row[0] == "":
                    i -= 1
                    continue
                if i == 0:
                    int_start = int(int_row[1]) + 1
                    int_end = int(int_row[2]) + 1
                    if int_start <= vcf_pos and vcf_pos < int_end:
                        f2.write("22" + "\t" + str(vcf_pos) + "\t" + "\n")
                        # print(vcf_pos)
                        # print(int_row[1] + "\t" + vcf_pos + "\t" + int_row[2])
                        i = BED_len - 1
                        break
                    if int_start <= vcf_pos and int_end <= vcf_pos:
                        # print("no" + "\t" + str(int_start) + "\t" + str(int_end))
                        i = BED_len - 1
                        break
                    else:
                        i = BED_len - 1
                        break
                else:
                    int_start = int(int_row[1]) + 1
                    int_end = int(int_row[2]) + 1
                    if int_start <= vcf_pos and vcf_pos < int_end:
                        print("yes" + "\t" + str(int_start) + "\t" + str(int_end))
                        f2.write("22" + "\t" + str(vcf_pos) + "\t" + "\n")
                        # print(vcf_pos)
                        # print(int_row[1] + "\t" + vcf_pos + "\t" + int_row[2])
                        i = BED_len - 1
                        break
                    if int_start <= vcf_pos and int_end <= vcf_pos:
                        # print("no" + "\t" + str(int_start) + "\t" + str(int_end))
                        i = BED_len - 1
                        break
                    else:
                        i -= 1
                        continue