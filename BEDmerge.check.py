# this script checks if all intervals were merged and there's no overlapping intervals in a BED file.

fn1 = "Chr22.regulation.cCREs.merged1.txt"
f1 = open(fn1, "r")

intervals = f1.readlines()
num_line = len(intervals)
# print(intervals)
# print(num_line)
i = 0

while True:
  if i >= num_line - 1:
    break
  else:
    line1 = intervals[i].split("\t")
    line2 = intervals[i+1].split("\t")
    # print(i)
    if line1[0] != "chr22":
      i +=1
      continue
    if line1[0] == "chr22":
      start2 = int(line2[1])
      end1 = int(line1[2])
      if end1 < start2:
        # print(str(i+1) + "\t" + str(start1))
        i +=1
        continue
      else:
        print("overlapping" + "\t" + str(i+1))
        i +=1
        continue