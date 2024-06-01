# this script makes a file containing regulatory region intervals in the reversed order from the original bed files

from file_read_backwards import FileReadBackwards

fn1 = "Chr22.regulation.cCREs.merged1.txt"
f1_in = FileReadBackwards(fn1, encoding="utf-8")

fn2 = "Chr22.regulation.cCREs.reversed.txt"
f2_in = open(fn2, "w")

while True:
  rev_row = f1_in.readline()
  print(rev_row)
  if rev_row == "":
    break
  else:
    if rev_row[0] == "c":
      f2_in.write(rev_row)
    else:
      break