# this script combines the lists of regulatory regions into one file. 

fn1 = "1kG.chr22.regulation.list.cCREs.txt"
f1 = open(fn1, "r")

fn2 = "1kG.chr22.regulation.list.ORegAnno.txt"
f2 = open(fn2, "r")

fn3 = "1kG.chr22.regulation.list.RefSeq.txt"
f3 = open(fn3, "r")

fn4 = "1kG.chr22.regulation.list.txt"
f4 = open(fn4, "w")

final_list = []

def pos_extraction(filename):
  while True:
    row = filename.readline().split(sep="\t")
    # print(row)
    if row == [""]:
      final_list.sort()
      break
    else:
      reg_pos = row[1]
      if reg_pos in final_list:
        continue
      if reg_pos not in final_list:
        final_list.append(reg_pos)
        continue

pos_extraction(f1)
pos_extraction(f2)
pos_extraction(f3)

# print(final_list)

for pos in final_list:
  f4.write("22" + "\t" + pos + "\t" + "\n")