# this script extracts all the types of annotation from snpEff annotation and consequence from VEP annotation.

import pysam

fn1 = "1kG-chr22-hg19-snpEff.vcf"
f1 = pysam.VariantFile(fn1)

fn2 = "1kG-chr22-hg19-vep.vcf"
f2 = pysam.VariantFile(fn2)

# fn2 = "variant_effect_output.txt"
# f2 = open(fn2, "r")

fn3 = "snpEff-VEP-types.txt"
f3 = open(fn3, "w")

conseq_list = []
anno_list = []

for record_snp in f1:
    try:
        snpEff = record_snp.info["ANN"]
        anno = snpEff[0].split(sep="|")[1]
        # print(anno)
        if anno in anno_list:
            continue
        if anno not in anno_list:
            anno_list.append(anno)
    except KeyError:
        continue

anno_list.sort()

for record_vep in f2:
    # print(record_vep)
    try:
        vep = record_vep.info["CSQ"]
        # print(vep)
        conseq = vep[0].split(sep="|")[1]
        # print(conseq)
        if conseq in conseq_list:
            continue
        if conseq not in conseq_list:
            conseq_list.append(conseq)
    except KeyError:
        # print(vep[0])
        continue

# while True:
#     line = f2.readline()
#     if len(line) == 0:
#         break
#     if line[0] == "#":
#         continue
#     else:
#         data = line.split(sep="\t")
#         conseq = data[6]
#         if conseq in conseq_list:
#             continue
#         if conseq not in conseq_list:
#             conseq_list.append(conseq)

conseq_list.sort()

if len(conseq_list) >= len(anno_list):
    f3.write("VEP_CONSEQ" + "\t" + "SNPEFF_ANN" + "\n")
    for i in range(len(conseq_list)):
        try:
            f3.write(conseq_list[i] + "\t" + anno_list[i] + "\n")
        except IndexError:
            f3.write(conseq_list[i] + "\n")

if len(conseq_list) < len(anno_list):
    f3.write("SNPEFF_ANN" + "\t" + "VEP_CONSEQ" + "\n")
    for i in range(len(anno_list)):
        try:
            f3.write(anno_list[i] + "\t" + conseq_list[i] + "\n")
        except IndexError:
            f3.write(anno_list[i] + "\n")