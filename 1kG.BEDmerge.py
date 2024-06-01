# this script sorts and merges the overlapping intervals into the one

fn1 = "Chr22.regulation.cCREs.txt"
f1 = open(fn1, "r")

fn2 = "Chr22.regulation.cCREs.merged1.txt"
f2 = open(fn2, "w")

intervals = f1.readlines()
num_line = len(intervals)
# print(intervals)
# print(num_line)
i = 0

while True:
  if i == num_line - 1:
    line = intervals[i].split("\t")
    start = int(line[1])
    end = int(line[2])
    f2.write("chr22" + "\t" + str(start) + "\t" + str(end) + "\t")
    break
  if i >= num_line:
    break
  else:
    line1 = intervals[i].split("\t")
    line2 = intervals[i+1].split("\t")
    # print(i)
    if line1[0] != "chr22":
      i +=1
      continue
    if line1[0] == "chr22":
      start1 = int(line1[1])
      start2 = int(line2[1])
      end1 = int(line1[2])
      end2 = int(line2[2])
      # print(line1[1] + line1[2] + "\t" + line2[1] + line2[2])
      if end1 < start2:
        f2.write("chr22" + "\t" + str(start1) + "\t" + str(end1) + "\t" + "\n")
        # print(str(i+1) + "\t" + str(start1))
        i +=1
        continue
      if end1 <= end2:
        f2.write("chr22" + "\t" + str(start1) + "\t" + str(end2) + "\t" + "\n")
        print(str(i+1) + "\t" + str(start1))
        i +=2
        continue
      if end2 <= end1:
        f2.write("chr22" + "\t" + str(start1) + "\t" + str(end1) + "\t" + "\n")
        print(str(i+1) + "\t" + str(start1))
        i +=2
        continue