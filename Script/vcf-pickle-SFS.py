# this script porduces the SFS file from pickle file that contains specific factors and its values.

import pickle
import numpy as np

def pickle_SFS(input_fn, output_fn, nc):
  input = pickle.load(open(input_fn, "rb"))
  output = open(output_fn, 'w')
  for factor in input:
    if "&" in factor:
      continue
    output.write(f"{factor}\n")
    sfs = np.zeros(nc+1, dtype=np.float32)
    for acan in input[factor]:
      ac = int(acan.split(sep="_")[1])
      sfs[ac] += input[factor][acan]
    sfs = sfs.astype(str)
    out_line = " ".join(sfs)
    output.write(out_line + "\n\n")

input_file = "Document/Pickle/1kG_sift4g_single_conseq.p"
output_file = "1kG-sift4g-single.txt"
nc = 5008

pickle_SFS(input_file, output_file, nc)