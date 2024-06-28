# this script sorts and merges the overlapping intervals into the one

fn1 = "Chr22-reg-ORegAnno.tsv"

fn2 = "Chr22-reg-ORegAnno-merged.tsv"

def BED_merging(input_fn, output_fn):
  input = open(input_fn, "r")
  intervals = input.readlines()
  num_line = len(intervals)
  input.close()
  output = open(output_fn, "w")
  i = 0 
  while True:
    if i == num_line - 1:
      line = intervals[i].split("\t")
      start = int(line[1])
      end = int(line[2])
      output.write("chr22" + "\t" + str(start) + "\t" + str(end) + "\t")
      output.close()
      break
    if i >= num_line:
      output.close()
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
          output.write("chr22" + "\t" + str(start1) + "\t" + str(end1) + "\t" + "\n")
          # print(str(i+1) + "\t" + str(start1))
          i +=1
          continue
        if end1 <= end2:
          output.write("chr22" + "\t" + str(start1) + "\t" + str(end2) + "\t" + "\n")
          print(str(i+1) + "\t" + str(start1))
          i +=2
          continue
        if end2 <= end1:
          output.write("chr22" + "\t" + str(start1) + "\t" + str(end1) + "\t" + "\n")
          print(str(i+1) + "\t" + str(start1))
          i +=2
          continue

def BED_merge_check(input_fn):
  input = open(input_fn, "r")
  intervals = input.readlines()
  num_line = len(intervals)
  overlap_list = []
  # print(intervals)
  # print(num_line)
  i = 0
  while True:
    if i >= num_line - 1:
      input.close()
      return overlap_list
      break
    if i >= num_line:
      input.close()
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
          overlap_list.append(str(i+1))
          i +=1
          continue


while True:
  BED_merging(fn1, fn2)
  overlap_list = BED_merge_check(fn2)
  if len(overlap_list) == 0:
    print("completed")
    break
  if len(overlap_list) != 0:
    print("on progress")
    fn1 = fn2
    continue