# this script sorts and merges the overlapping intervals into the one
# this script relies on the intervals order. startpoints have to be sorted in increasing order.

fn1 = "Chr22.exon.RefSeq.sorted12.tsv"
fn2 = "Chr22.exon.RefSeq.merged12.tsv"

def BED_merging(input_fn, output_fn):
  input = open(input_fn, "r")
  input_row = input.readlines()
  num_line = len(input_row)
  input.close()
  output = open(output_fn, "w")
  i = 0 
  sample_line = input_row[1].replace("\n", "").split(sep="\t")
  if len(sample_line) >= 6:
    while True:
      if i == num_line - 1:
        line = input_row[i].replace("\n", "").split("\t")
        start = int(line[1])
        end = int(line[2])
        output.write("chr22" + "\t" + str(start) + "\t" + str(end) + "\t" + line[5] + "\t" + "\n")
        i = 0
        output.close()
        break
      if i >= num_line - 1:
        output.close()
        i = 0
        break
      else:
        line1 = input_row[i].replace("\n", "").split("\t")
        line2 = input_row[i+1].replace("\n", "").split("\t")
        # print(line1)
        if line1[0] != "chr22":
          i +=1
          continue
        if line1[0] == "chr22":
          start1 = int(line1[1])
          start2 = int(line2[1])
          end1 = int(line1[2])
          end2 = int(line2[2])
          if end1 < start2:
            # no merging
            output.write("chr22" + "\t" + str(start1) + "\t" + str(end1) + "\t" + line1[5] + "\t" + "\n")
            i +=1
            continue
          if end1 < end2:
            # merging in case of bigger second interval
            info1 = line1[5]
            info2 = line2[5]
            if info1 == info2:
              output.write("chr22" + "\t" + str(start1) + "\t" + str(end2) + "\t" + info1 + "\t" + "\n")
              i +=2
              continue
            if start1 < start2:
              output.write("chr22" + "\t" + str(start1) + "\t" + str(start2 - 1) + "\t" + info1 + "\t" + "\n")
              output.write("chr22" + "\t" + str(start2) + "\t" + str(end1) + "\t" + "+,-" + "\t" + "\n")
              output.write("chr22" + "\t" + str(end1 + 1) + "\t" + str(end2) + "\t" + info2 + "\t" + "\n")
              i +=2
              continue
            if start1 == start2:
              output.write("chr22" + "\t" + str(start2) + "\t" + str(end1) + "\t" + "+,-" + "\t" + "\n")
              output.write("chr22" + "\t" + str(end1 + 1) + "\t" + str(end2) + "\t" + info2 + "\t" + "\n")
              i +=2
              continue
          if end2 < end1:
            # merging in case of bigger first interval
            info1 = line1[5]
            info2 = line2[5]
            if info1 == info2:
              output.write("chr22" + "\t" + str(start1) + "\t" + str(end1) + "\t" + info1 + "\t" + "\n")
              i +=2
              continue
            if start1 < start2:
              output.write("chr22" + "\t" + str(start1) + "\t" + str(start2 - 1) + "\t" + info1 + "\t" + "\n")
              output.write("chr22" + "\t" + str(start2) + "\t" + str(end2) + "\t" + "+,-" + "\t" + "\n")
              output.write("chr22" + "\t" + str(end2 + 1) + "\t" + str(end1) + "\t" + info1 + "\t" + "\n")
              i +=2
              continue
            if start1 == start2:
              output.write("chr22" + "\t" + str(start2) + "\t" + str(end2) + "\t" + "+,-" + "\t" + "\n")
              output.write("chr22" + "\t" + str(end2 + 1) + "\t" + str(end1) + "\t" + info1 + "\t" + "\n")
              i +=2
              continue
          if end1 == end2:
            info1 = line1[5]
            info2 = line2[5]
            if start1 == start2:
              if info1 == info2:
                output.write("chr22" + "\t" + str(start1) + "\t" + str(end1) + "\t" + info1 + "\t" + "\n")
                i +=2
                continue
              else:
                output.write("chr22" + "\t" + str(start1) + "\t" + str(end1) + "\t" + "+,-" + "\t" + "\n")
                i +=2
                continue
            else:
              if info1 == info2:
                output.write("chr22" + "\t" + str(start1) + "\t" + str(end1) + "\t" + info1 + "\t" + "\n")
                i +=2
                continue
              else:
                output.write("chr22" + "\t" + str(start1) + "\t" + str(start2 - 1) + "\t" + info1 + "\t" + "\n")
                output.write("chr22" + "\t" + str(start2) + "\t" + str(end1) + "\t" + "+,-" + "\t" + "\n")
                i +=2
                continue
  else:
    while True:
      if i == num_line - 1:
        line = input_row[i].split("\t")
        start = int(line[1])
        end = int(line[2])
        output.write("chr22" + "\t" + str(start) + "\t" + str(end) + "\t" + line[3] + "\t" + "\n")
        i = 0
        output.close()
        break
      if i >= num_line - 1:
        output.close()
        i = 0
        break
      else:
          line1 = input_row[i].split("\t")
          line2 = input_row[i+1].split("\t")
          # print(line1)
          # print(line2)
          start1 = int(line1[1])
          start2 = int(line2[1])
          end1 = int(line1[2])
          end2 = int(line2[2])
          if end1 < start2:
            # no merging
            output.write("chr22" + "\t" + str(start1) + "\t" + str(end1) + "\t" + line1[3] + "\t" + "\n")
            # print(str(i+1) + "\t" + str(start1))
            i +=1
            continue
          if end1 < end2:
            # merging in case of bigger second interval
            info1 = line1[3]
            info2 = line2[3]
            if info1 == info2:
              output.write("chr22" + "\t" + str(start1) + "\t" + str(end2) + "\t" + info1 + "\t" + "\n")
              i +=2
              continue
            if start1 < start2:
              output.write("chr22" + "\t" + str(start1) + "\t" + str(start2 - 1) + "\t" + info1 + "\t" + "\n")
              output.write("chr22" + "\t" + str(start2) + "\t" + str(end1) + "\t" + "+,-" + "\t" + "\n")
              output.write("chr22" + "\t" + str(end1 + 1) + "\t" + str(end2) + "\t" + info2 + "\t" + "\n")
              i +=2
              continue
            if start1 == start2:
              output.write("chr22" + "\t" + str(start2) + "\t" + str(end1) + "\t" + "+,-" + "\t" + "\n")
              output.write("chr22" + "\t" + str(end1 + 1) + "\t" + str(end2) + "\t" + info2 + "\t" + "\n")
              i +=2
              continue
          if end2 < end1:
            # merging in case of bigger first interval
            info1 = line1[3]
            info2 = line2[3]
            if info1 == info2:
              output.write("chr22" + "\t" + str(start1) + "\t" + str(end1) + "\t" + info1 + "\t" + "\n")
              i +=2
              continue
            if start1 < start2:
              output.write("chr22" + "\t" + str(start1) + "\t" + str(start2 - 1) + "\t" + info1 + "\t" + "\n")
              output.write("chr22" + "\t" + str(start2) + "\t" + str(end2) + "\t" + "+,-" + "\t" + "\n")
              output.write("chr22" + "\t" + str(end2 + 1) + "\t" + str(end1) + "\t" + info1 + "\t" + "\n")
              i +=2
              continue
            if start1 == start2:
              output.write("chr22" + "\t" + str(start2) + "\t" + str(end2) + "\t" + "+,-" + "\t" + "\n")
              output.write("chr22" + "\t" + str(end2 + 1) + "\t" + str(end1) + "\t" + info1 + "\t" + "\n")
              i +=2
              continue
          if end1 == end2:
            info1 = line1[3]
            info2 = line2[3]
            if start1 == start2:
              if info1 == info2:
                output.write("chr22" + "\t" + str(start1) + "\t" + str(end1) + "\t" + info1 + "\t" + "\n")
                i +=2
                continue
              else:
                output.write("chr22" + "\t" + str(start1) + "\t" + str(end1) + "\t" + "+,-" + "\t" + "\n")
                i +=2
                continue
            else:
              if info1 == info2:
                output.write("chr22" + "\t" + str(start1) + "\t" + str(end1) + "\t" + info1 + "\t" + "\n")
                i +=2
                continue
              else:
                output.write("chr22" + "\t" + str(start1) + "\t" + str(start2 - 1) + "\t" + info1 + "\t" + "\n")
                output.write("chr22" + "\t" + str(start2) + "\t" + str(end1) + "\t" + "+,-" + "\t" + "\n")
                i +=2
                continue

def BED_merge_check(input_fn):
  input = open(input_fn, "r")
  intervals = input.readlines()
  num_line = len(intervals)
  overlap_list = []
  # print(intervals)
  # print(num_line)
  g = 0
  while True:
    if g >= num_line - 1:
      input.close()
      return overlap_list
    else:
      line1 = intervals[g].split("\t")
      line2 = intervals[g+1].split("\t")
      if len(line2) == 0:
        input.close()
        return overlap_list
      else:
        # print(line2)
        # print(i)
        if line1[0] != "chr22":
          g +=1
          continue
        if line1[0] == "chr22":
          start2 = int(line2[1])
          end1 = int(line1[2])
          if end1 < start2:
            # print(str(i+1) + "\t" + str(start1))
            g +=1
            continue
          else:
            # print("overlapping" + "\t" + str(g+1))
            overlap_list.append(str(g+1))
            g +=1
            continue


out_list = BED_merge_check(fn1)
if len(out_list) == 0:
  print("completed")
if len(out_list) != 0:
  print("on progress")
  BED_merging(fn1, fn2)
  fn1 = fn2
  # startpoint(fn2, fn2)