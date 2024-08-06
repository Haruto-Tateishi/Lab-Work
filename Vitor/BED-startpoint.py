# this script checks if the starting points of the intervals of BED files are sorted or not. If not, sort it

fn1 = "Chr22-CpG-merged.tsv"
fn2 = "a.tsv"


def start_sort(list_in, len_in):
  for i in range(len_in):
  # for i in range(1000):
    BED_row = list_in[i].split(sep="\t")
    if BED_row[0] != "chr22":
      continue
    else:
      start = int(BED_row[1])
      current_point = i
      for g in range(1, len_in):
        try:
          comp_row = list_in[i + g].split(sep="\t")
          comp_start = int(comp_row[1])
          if start > comp_start:
            list_in[current_point] = "\t".join(comp_row)
            list_in[i+g] = "\t".join(BED_row)
            print(str(start) + "\t" + str(comp_start))
            current_point = i + g
        except IndexError:
          break

def start_check(count, list_in, len_in):
  for t in range(len_in - 1):
    line1 = list_in[t].split("\t")
    line2 = list_in[t+1].split("\t")
    if line1[0] == "chr22":
      start1 = int(line1[1])
      start2 = int(line2[1])
      # if start1 < start2:
      #   print("yes")
      if start1 > start2:
        # print("no")
        # print(str(t) + "\t" + str(start1) + "\t" + str(start2))
        count +=1
  return count


def startpoint(in_file, out_file):
  f_in = open(in_file, "r")
  BED_lines = f_in.readlines()
  BED_len = len(BED_lines)
  f_in.close()
  while True:
    count_in = 0
    count_out = start_check(count_in, BED_lines, BED_len)
    if count_out == 0:
      print("completed")
      break
    if count_out > 0:
      print("on progress")
      # print(count_in)
      # print(count_out)
      start_sort(BED_lines, BED_len)
      continue

  f_out = open(out_file, "w")
  for h in BED_lines:
    f_out.write(h)
  f_out.close()


startpoint(fn1, fn2)