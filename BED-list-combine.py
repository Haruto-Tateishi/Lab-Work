# this script combines the lists of regulatory regions into one file. 

fn1 = "Chr22-reg-list-ENCODE.tsv"
f1 = open(fn1, "r")

fn2 = "Chr22-reg-list-ORegAnno.tsv"
f2 = open(fn2, "r")

fn3 = "Chr22-reg-list-Hancer.tsv"
f3 = open(fn3, "r")

fn4 = "Chr22-reg-list.tsv"
f4 = open(fn4, "w")

pos_dic = {}

def pos_extraction(filename, abbreviation):
  while True:
    row = filename.readline().split(sep="\t")
    print(row)
    if row == [""]:
      break
    else:
      reg_pos = row[1]
      if reg_pos in pos_dic.keys():
        pos_sup = pos_dic[reg_pos] + "," + abbreviation
        pos_dic[reg_pos] = pos_sup
        continue
      if reg_pos not in pos_dic.keys():
        pos_dic[reg_pos] = abbreviation
        continue

pos_extraction(f1, "EN")
pos_extraction(f2, "OR")
pos_extraction(f3, "HA")

final_dic = dict(sorted(pos_dic.items()))

# print(final_dic)

for key in final_dic.keys():
  f4.write("22" + "\t" + key + "\t" + final_dic[key] + "\t" + "\n")