# deal with vcf file by using bed file.
# For example, we can use bed file to filter vcf file.
# We can also use bed file to get the position of vcf file.

import pickle
import os 
import logging

# list all files in a directory
def list_documents(directory, extension):
    try:
        file_names = os.listdir(directory)
        file_names = sorted(file_names)
        document_extensions = {extension}
        documents = [f for f in file_names if os.path.isfile(os.path.join(directory, f)) and os.path.splitext(f)[1] in document_extensions]
        return documents
    except Exception as e:
        logging.error(f"Error listing documents: {e}", exc_info=True)

# read bed file and fetch lines
def bed_file_fetch(bed_fn, chr):
    with open(bed_fn, 'r') as f:
        lines = f.readlines()
    new_lines = []
    for line in lines:
        cmp_chr = line.strip().split(sep="\t")[0]
        if cmp_chr == chr:
            new_lines.append(line)
        else:
            continue
    return new_lines

# combine all pos lists
def pos_list_generate(input_fn, chr, conseq, pos_list):
  with open(input_fn, "rb") as f:
    ref_dict = pickle.load(f)
  for GC_bias in ref_dict.keys():
      if conseq == "syn":
        for codon in ref_dict[GC_bias].keys():
            for bin_num in ref_dict[GC_bias][codon].keys():
              cur_pos_list = list(ref_dict[GC_bias][codon][bin_num].keys())
              pos_list = pos_list + cur_pos_list
              pass
      elif conseq == "int":
          for bin_num in ref_dict[GC_bias].keys():
              cur_pos_list = list(ref_dict[GC_bias][bin_num].keys())
              pos_list = pos_list + cur_pos_list
              pass
#   print(pos_list)
  return sorted(pos_list)

# make super dictinary
# chr
#   pos
#     bed interval
def super_dict_generate(pos_list, bed_lines, super_dict, chr):
  zero_one = [0,1]
  super_dict[chr] = {}
  for pos in pos_list:
      pos = int(pos)
      bed_interval = False
      for bed_line in bed_lines:
          start_pos = int(bed_line.strip().split(sep="\t")[1]) + 1
          end_pos = int(bed_line.strip().split(sep="\t")[2]) + 1
        #   print(start_pos, end_pos, pos)
          if pos >= start_pos and pos <= end_pos:
              bed_interval = True
              break
          elif pos < start_pos:
              print(pos)
              break            
          else:
              bed_lines.remove(bed_line)
              continue
      super_dict[chr][pos] = bed_interval
      pass
  return super_dict

# main function
# 1. list all files in the directory
# 2. read bed file and fetch lines
# 3. combine all pos lists
# 4. make super dictinary
# 5. save super dictinary
bank = "1kG"
directory = "Document/Pickle/MutationRate/SuperDictionary/1kG_GC_syn_int"
documents = list_documents(directory, extension=".p")
bed_dir = "Document/Bedfile/MutationRate/mask_strict.bed"
super_dict = {}
for chr in range(1, 23):
  print(f"chr{chr} started")
  std_chr = f"chr{chr}"
  bed_lines = bed_file_fetch(bed_dir, std_chr)
  pos_list = []
  for input_fn in documents:
      inp_chr = input_fn.split(sep="_")[1]
      if inp_chr == std_chr:
          conseq = input_fn.split(sep="_")[3].split(sep=".")[0]
          input_fn = os.path.join(directory, input_fn)
          pos_list = pos_list_generate(input_fn, inp_chr, conseq, pos_list)
          pass
      else:
          continue
  super_dict = super_dict_generate(pos_list, bed_lines, super_dict, std_chr)
  print("end")
pickle_fn = f"{bank}_super_dict_bed_strict.p"
with open(pickle_fn, 'wb') as f:
    pickle.dump(super_dict, f)