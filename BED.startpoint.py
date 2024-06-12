# this script checks if the starting points of the intervals of BED files are sorted or not

fn1 = "Chr22.CpG.tsv"
f1 = open(fn1, "r")


intervals = f1.readlines()
num_line = len(intervals)


for i in range(num_line - 1):
  line1 = intervals[i].split("\t")
  line2 = intervals[i+1].split("\t")
  if line1[0] == "chr22":
    start1 = int(line1[1])
    start2 = int(line2[1])
    # if start1 < start2:
    #   print("yes")
    if start1 > start2:
      # print("no")
      print(str(i) + "\t" + str(start1) + "\t" + str(start2))