# This script generates SFS txt files by analyzing the pickle file that contains the dictionary of flanking base, consequence, location and allele count of each consequences.

import pickle
import os.path as op
import numpy as np

# Make a list of all possible CpG combinations
def CpG_list_function():
  CpG_list = []
  codon_list = ["CG/A","CG/C","CG/T","C/AG","C/GG","C/GT","CA/G","CC/G","CT/G","A/CG","G/CG","T/CG","GC/A","GC/G","GC/T","G/AC","G/CC","G/TC","GA/C","GG/C","GT/C","A/GC","C/GC","T/GC"]
  base_list = ["A","C","G","T"]
  for codon in codon_list:
    index = codon.find("/")
    if index == 2:
      for base in base_list:
        CpG = codon + base
        if CpG in CpG_list:
          continue
        else:
          CpG_list.append(CpG)
          pass
    elif index == 1:
      for base in base_list:
        CpG = base + codon
        if CpG in CpG_list:
          continue
        else:
          CpG_list.append(CpG)
          pass
  return CpG_list

# Function to generate SFS from the pickle file
# The function takes the pickle file name, output file name, control and experimental consequence names, and the number of alleles as input
# and generates SFS for each flanking base
# The function also filters out CpG sites
def flanking_pickle_SFS_single(pickle_fn, output_fn, con, exp, nc):
  keyerror_skipped = 0
  super_dict = pickle.load(open(pickle_fn, 'rb'))
  CpG_list = CpG_list_function()
  result_dict ={} # flanking base as first key, controll/experimental as second key, and list of SFS as value.
  for chr in super_dict:
    for flanking in super_dict[chr]:
      if flanking in CpG_list:
        continue
      elif flanking not in result_dict:
        result_dict[flanking] = {}
        result_dict[flanking][con] = np.zeros(nc+1, dtype=np.float32)
        result_dict[flanking][exp] = np.zeros(nc+1, dtype=np.float32)
      try:
        con_dict = super_dict[chr][flanking][con]
        con_loc_list = list(con_dict.keys())
        exp_dict = super_dict[chr][flanking][exp]
      except KeyError:
        keyerror_skipped+=1
        continue
      for exp_loc in exp_dict:
        exp_ac = exp_dict[exp_loc]
        con_loc = min(con_loc_list, key=lambda x: abs(x - exp_loc))
        con_loc_list.remove(con_loc)
        con_ac = con_dict[con_loc]
        result_dict[flanking][con][con_ac] +=1
        result_dict[flanking][exp][exp_ac] +=1
  print(f"SFS generation is done. Skippeed {keyerror_skipped} times")
  # output the result dictionary
  output = open(output_fn, "w")
  for flanking in result_dict:
    output.write(f"{flanking}\n")
    con_sfs = result_dict[flanking][con]
    con_sfs = con_sfs.astype(str)
    con_sfs = " ".join(con_sfs)
    exp_sfs = result_dict[flanking][exp]
    exp_sfs = exp_sfs.astype(str)
    exp_sfs = " ".join(exp_sfs)
    output.write(f"{con}\n")
    output.write(f"{con_sfs}\n\n")
    output.write(f"{exp}\n")
    output.write(f"{exp_sfs}\n\n")
  output.close()
  print("Done")

# Function to generate SFS from the pickle file
# The function takes the pickle file name, output file name, control and experimental consequence names, and the number of alleles as input
# and generates SFS for each flanking base
# The function also filters out CpG sites
# The function is similar to flanking_pickle_SFS_single, but it generates SFS for all flanking bases
def flanking_pickle_SFS_all(pickle_fn, output_fn, con, exp, nc):
  keyerror_skipped = 0
  super_dict = pickle.load(open(pickle_fn, 'rb'))
  # CpG_list = CpG_list_function()
  result_dict ={} # controll/experimental as first key, and list of SFS as value.
  result_dict[con] = np.zeros(nc+1, dtype=np.float32)
  result_dict[exp] = np.zeros(nc+1, dtype=np.float32)
  for chr in super_dict:
    for flanking in super_dict[chr]:
      # if flanking in CpG_list:
      #   continue
      try:
        con_dict = super_dict[chr][flanking][con]
        con_loc_list = list(con_dict.keys())
        exp_dict = super_dict[chr][flanking][exp]
      except KeyError:
        keyerror_skipped+=1
        continue
      for exp_loc in exp_dict:
        exp_ac = exp_dict[exp_loc]
        con_loc = min(con_loc_list, key=lambda x: abs(x - exp_loc))
        con_loc_list.remove(con_loc)
        con_ac = con_dict[con_loc]
        result_dict[con][con_ac] +=1
        result_dict[exp][exp_ac] +=1
  print(f"SFS generation is done. Skippeed {keyerror_skipped} times")
  # output the result dictionary
  output = open(output_fn, "w")
  con_sfs = result_dict[con]
  con_sfs = con_sfs.astype(str)
  con_sfs = " ".join(con_sfs)
  exp_sfs = result_dict[exp]
  exp_sfs = exp_sfs.astype(str)
  exp_sfs = " ".join(exp_sfs)
  output.write(f"{con}\n")
  output.write(f"{con_sfs}\n\n")
  output.write(f"{exp}\n")
  output.write(f"{exp_sfs}\n\n")
  output.close()
  print("Done")

# Function to generate SFS from the pickle file
# The function takes the pickle file name, output file name, control and experimental consequence names, and the number of alleles as input
# and generates SFS for each flanking base
# The function also filters out CpG sites
# The function is similar to flanking_pickle_SFS_single, but it generates SFS for all flanking bases
# it saves the SFS in pickle file
def flanking_pickle_all_pickle(pickle_fn, output_fn, con, exp, nc):
  keyerror_skipped = 0
  super_dict = pickle.load(open(pickle_fn, 'rb'))
  CpG_list = CpG_list_function()
  result_dict ={} # controll/experimental as first key, and AC_#_AN_# as second key and counts as value.
  result_dict[con] = {}
  result_dict[exp] = {}
  for chr in super_dict:
    for flanking in super_dict[chr]:
      if flanking in CpG_list:
        continue
      try:
        con_dict = super_dict[chr][flanking][con]
        con_loc_list = list(con_dict.keys())
        exp_dict = super_dict[chr][flanking][exp]
      except KeyError:
        keyerror_skipped+=1
        continue
      for exp_loc in exp_dict:
        exp_ac = exp_dict[exp_loc]
        exp_acan = "AC_" + str(exp_ac) + "_AN_7562"
        con_loc = min(con_loc_list, key=lambda x: abs(x - exp_loc))
        con_loc_list.remove(con_loc)
        con_ac = con_dict[con_loc]
        con_acan = "AC_" + str(con_ac) + "_AN_7562"
        if exp_acan in result_dict[exp].keys():
          result_dict[exp][exp_acan] += 1
          pass
        else:
          result_dict[exp][exp_acan] = 1
          pass
        if con_acan in result_dict[con].keys():
          result_dict[con][con_acan] += 1
          pass
        else:
          result_dict[con][con_acan] = 1
          pass
  print(f"Done. Skippeed {keyerror_skipped} times")
  # output the result dictionary
  pickle.dump(result_dict, open(output_fn, 'wb'))
  print("Done")

# Function to generate SFS from the pickle file
# it uses synonymnous codon pair as key of super dictionary
def flanking_syn_pair_pickle_SFS_all(pickle_fn, output_fn, con, exp, nc):
  processed = 0
  keyerror_skipped = 0
  CpG_skipped = 0
  super_dict = pickle.load(open(pickle_fn, 'rb'))
  CpG_list = CpG_list_function()
  result_dict ={} # controll/experimental as first key, and AC_#_AN_# as second key and counts as value.
  for chr in super_dict:
    for flanking in super_dict[chr]:
      if flanking in CpG_list:
        CpG_skipped +=1
        continue
      try:
        con_dict = super_dict[chr][flanking][con]
        con_loc_list = list(con_dict.keys())
        exp_dict = super_dict[chr][flanking][exp]
      except KeyError:
        keyerror_skipped+=1
        continue
      for syn_codon in exp_dict.keys():
        for exp_loc in exp_dict[syn_codon]:
          exp_ac = exp_dict[syn_codon][exp_loc]
          exp_acan = "AC_" + str(exp_ac) + "_AN_" + str(nc)
          con_loc = min(con_loc_list, key=lambda x: abs(x - exp_loc))
          con_loc_list.remove(con_loc)
          con_ac = con_dict[con_loc]
          con_acan = "AC_" + str(con_ac) + "_AN_" + str(nc)
          processed +=1
          if syn_codon in result_dict.keys():
            if exp_acan in result_dict[syn_codon][exp].keys():
              result_dict[syn_codon][exp][exp_acan] += 1
              pass
            else:
              result_dict[syn_codon][exp][exp_acan] = 1
              pass
            if con_acan in result_dict[syn_codon][con].keys():
              result_dict[syn_codon][con][con_acan] += 1
              pass
            else:
              result_dict[syn_codon][con][con_acan] = 1
              pass
          elif syn_codon not in result_dict.keys():
            result_dict[syn_codon] = {}
            result_dict[syn_codon][con] = {}
            result_dict[syn_codon][exp] = {}
            result_dict[syn_codon][exp][exp_acan] = 1
            result_dict[syn_codon][con][con_acan] = 1
            pass
  print(f"Done. Processed {processed} time. Skippeed key_error:{keyerror_skipped} times, CpG:{CpG_skipped} times")
  # output the result dictionary
  pickle.dump(result_dict, open(output_fn, 'wb'))
  print("Done")


# Example usage
pickle_fn = "Document/Pickle/FlankingBases/syn_codon/1kG_vep_pickorder_flanking_syn_codon.p"
output_fn = "Document/Pickle/1kG_vep_pickorder_flanking_CpGfiltered_syn_codon.p"
nc = 5008
flanking_syn_pair_pickle_SFS_all(pickle_fn, output_fn, "intergenic_variant", "synonymous_variant", nc)