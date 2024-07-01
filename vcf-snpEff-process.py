# this script generates the text file containing not only number of SNPs but also annotations

fn1 = "chr22.annotation.vcf"
f1_in = open(fn1, 'rb')

fn2 = "integrated_call_samples_v3.20130502.ALL.panel.txt"
f2_in = open(fn2, 'r')

fn3 = "1kG.chr22.annotation.vcf"
f3_in = open(fn3, 'w')


lines = f2_in.readlines()
id_dic = {}
pop_dic = {}
for line in lines[1:]:
    temp = line.split()
    id_dic[temp[0]] = temp[1]
    if temp[1] not in pop_dic:
        pop_dic[temp[1]] = [temp[0]]
    else:
        pop_dic[temp[1]].append(temp[0])
pop_list = pop_dic.keys()

# print(id_dic)
# print(pop_dic)
f2_in.close

# print(pop_list)

while True:
# for g in range(500):
    line = f1_in.readline().decode()
    if len(line) <= 1:
      break
    if line[1] == "#":
        f3_in.write(line)
    if line[1] == "C":
       sample_id_list = line.strip().split(sep="\t")[9:]
       top_list = []
       for i1 in range(9):
           top_list.append(line.strip().split(sep="\t")[i1])
       for i2 in pop_list:
           vari_pop = "VARI_" + i2
           top_list.append(vari_pop)
       top_str = "\t".join(top_list)
       f3_in.write(top_str + "\n")
      #  print(top_str)
    if line[0] != "#":
      row = line.strip().split(sep="\t")
      ref_alt = row[3:5]
      len_ref_alt = 0
      for i3 in ref_alt:
          len_ref_alt = len_ref_alt + len(i3)
      if len_ref_alt ==2:
          row_data = row[9:]
          sample_id = []
          sample_pop = []
          num_sample = len(row_data)
          for i4 in range(num_sample):
            if row_data[i4][0] =='1' or row_data[i4][2] =='1':
                sample_id.append(sample_id_list[i4])
            if row_data[i4][0] =='1' and row_data[i4][2] =='1':
                sample_id.append(sample_id_list[i4])
          # print(sample_id)
          for i5 in sample_id:
              sample_pop.append(id_dic.get(i5))
          # print(sample_pop)
          data_list = []
          for i6 in range(9):
              data_list.append(row[i6])
          for i7 in pop_list:
              data_list.append(str(sample_pop.count(i7)))
          data_str = "\t".join(data_list)
          f3_in.write(data_str + "\n")
          # print(data_str)

f3_in.close