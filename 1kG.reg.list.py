# this script makes a list of regulatory region SNPs based on BED file


fn1 = "1kG.chr22.annotation.txt"
f1_in = open(fn1, 'r')

fn2 = "1kG.chr22.regulation.list.cCREs.txt"
f2_in = open(fn2, 'w')

fni = "Chr22.regulation.cCREs.reversed.txt"
fni_in = open(fni, "r")


while True:
# for g in range(1000):
    vcf_data = f1_in.readline().split(sep="\t")
    if len(vcf_data) == 0:
        break 
    if vcf_data[0] == "22":
        vcf_pos = int(vcf_data[1])
        print(vcf_pos)
        if vcf_pos >= 50800000:
            break
        else:
            while True:
                int_row = fni_in.readline().split(sep="\t")
                # print(int_row)
                if int_row[0] == "":
                    fni_in.seek(0)
                    break   
                else:
                    int_start = int(int_row[1])
                    int_end = int(int_row[2])
                    if int_start <= vcf_pos and int_end >= vcf_pos:
                        print("yes")
                        f2_in.write("22" + "\t" + str(vcf_pos) + "\t" + "\n")
                        # print(vcf_pos)
                        # print(int_row[1] + "\t" + vcf_pos + "\t" + int_row[2])
                        fni_in.seek(0)
                        break
                    if int_start <= vcf_pos and int_end < vcf_pos:
                        print("no" + "\t" + str(int_start) + "\t" + str(int_end))
                        fni_in.seek(0)
                        break
                    