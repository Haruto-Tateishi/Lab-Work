# this script creates the SFS for each category of annotation/consequences.

import matplotlib.pyplot as plt
import pysam

fn1 = "1kG-chr22-hg19-vep.vcf"
# f1 = pysam.VariantFile(fn1)

fn2 = "VEP-types.txt"
f2 = open(fn2, "r")

# plot scatter graph for whole chromsome
# whole_list = []
# for line in f1:
#     allelec = line.info["AC"][0]
#     whole_list.append(allelec)
#     continue

# sorted_whole = sorted(whole_list)
# max_c = max(sorted_whole)
# list_new = []
# for j in range(1, max_c+1):
#     counted = sorted_whole.count(j)
#     list_new.append(counted)
# list_length = len(list_new)
# x1 = []
# y1 = []
# for k in range(1, list_length+1):
#     x1.append(k)
#     y_values = list_new[k-1]
#     y1.append(y_values)
#     continue
# plt.plot(x1, y1)
# plt.xlabel('allele counts')
# plt.ylabel('counts of position')
# plt.title('SFS of whole chrom')
# plt.show()
# f1.close()

f3 = pysam.VariantFile(fn1)

conseq_list = f2.readlines()
conseq_len = len(conseq_list)
bin_dict = {}
i = 1
while True:
    conseq_type = conseq_list[i].replace("\n", "").split(sep="\t")[0]
    bin_dict[conseq_type] = []
    i +=1
    if i == conseq_len:
        break
    else:
        continue

# import time

for record in f3:
    try:
        veps = record.info["CSQ"]
        anno_list = []
        for v in veps:
            conseq = v.split(sep="|")[1]
            anno_list.append(conseq)
        vep = sorted(set(anno_list))
        # print(vep)
        ac = record.info["AC"][0]
        # print(ac)
        # time.sleep(1/20)
        for annotation in vep:
            if "splice_donor_variant" in annotation:
                print(annotation, ac)
            bin_dict[annotation].append(ac)
    except KeyError:
        # print(annotation)
        continue

# print(bin_dict)

# # plot bar chart
# for conseq_a in bin_dict.keys():
#     sorted_list = sorted(bin_dict[conseq_a])
#     range = (1, 5008)
#     bins = 5008
#     plt.hist(sorted_list, bins, range, color = 'blue',
#         histtype = 'bar', rwidth = 0.8)
#     plt.xlabel('allele counts')
#     plt.ylabel('counts of position')
#     plt.title('SFS of ' + conseq_a)
#     plt.show()


# plot scatter graph
for conseq_a in bin_dict.keys():
    sorted_list = sorted(bin_dict[conseq_a])
    print(conseq_a)
    # print(sorted_list)
    max_count = max(sorted_list)
    new_list = []
    for g in range(1, max_count+1):
        counts = sorted_list.count(g)
        new_list.append(counts)
    bin_dict[conseq_a] = new_list


for conseq_b in bin_dict.keys():
    bin_list = bin_dict[conseq_b]
    list_len = len(bin_list)
    x = []
    y = []
    for h in range(1, list_len+1):
        x.append(h)
        y_value = bin_list[h-1]
        y.append(y_value)
        continue
    plt.plot(x, y)
    plt.xlabel('allele counts')
    plt.ylabel('counts of position')
    plt.title('SFS of ' + conseq_b)
    plt.show()


f3.close()
f2.close()