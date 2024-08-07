# This script combines the tables of AUC in this case 6 tables

super_list = [[] for _ in range(135)]

ref_file = "synonymous_codon_pairs.txt"
ref = open(ref_file, 'r')
ref_dict = {}
for g in range(67):
  ref_line = ref.readline()
  ref_conseq_1 = ref_line.split(sep="\t")[0]
  ref_conseq_2 = ref_line.split(sep="\t")[1]
  if ref_conseq_1 not in ref_dict:
    ref_aa = ref_line.strip().split(sep="\t")[2]
    ref_dict[ref_conseq_1] = ref_aa
  if ref_conseq_2 not in ref_dict:
    ref_aa = ref_line.strip().split(sep="\t")[2]
    ref_dict[ref_conseq_2] = ref_aa

fn1 = "1kG-SFS-syn-single-integrate-table.tsv"
f1 = open(fn1, 'w')

num_list = [30, 50, 100, 200, 500]
for nc in num_list:
  fn2 = f"1kG-SFS-downsample-integrate-syn-single-{nc}.tsv"
  f2 = open(fn2, 'r')
  lines = f2.readlines()
  super_list[0].append(f"{str(nc)}-conseq")
  super_list[0].append(f"{str(nc)}-aa")
  for h in range(1, 135):
    line = lines[h]
    codon_pair = line.split(sep="\t")[0]
    conseq = codon_pair[0:3]
    aa = ref_dict[conseq]
    super_list[h].append(codon_pair)
    super_list[h].append(aa)

fn2 = "1kG-SFS-syn-single-integrate.tsv"
f2 = open(fn2, 'r')
lines = f2.readlines()
super_list[0].append("5008-conseq")
super_list[0].append("5008-aa")
for j in range(1, 135):
    line = lines[j]
    codon_pair = line.split(sep="\t")[0]
    conseq = codon_pair[0:3]
    aa = ref_dict[conseq]
    super_list[j].append(codon_pair)
    super_list[j].append(aa)

for k in range(135):
  out_line = super_list[k]
  out_line.append("\n")
  out_line = "\t".join(out_line)
  f1.write(out_line)