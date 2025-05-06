# This script goes through the outputs of the SF_Ratios_run.py

import os
import numpy as np

# Function to list all documents in a directory
def list_documents(directory):
    try:
        # Get the list of all files in the directory
        file_names = os.listdir(directory)
        file_names = sorted(file_names)
        # Filter to include only documents (e.g., .txt, .pdf, .docx)
        document_extensions = {'.out'}  # Add more extensions if needed
        documents = [f for f in file_names if os.path.isfile(os.path.join(directory, f)) and os.path.splitext(f)[1] in document_extensions]

        return documents
    except Exception as e:
        print(f"Error: {e}")

# Function to prepare the output file
def output_prep(output_fn):
   output = open(output_fn, 'w')
   header = ["codon", "2Ns", "clm_ratio", '#_1st_bin<10', 'total_SNPs(w/o 0 bin)', "rev_codon", "rev_2Ns", "rev_clm_ratio", 'rev_#_1st_bin<10', 'rev_total_SNPs(w/o 0 bin)', 'Comparison', "codon", "2Ns", "rev_codon", "rev_2Ns", "closeness"]
   header = "\t".join(header) + "\n"
   output.write(header)
   output.close()
  
# Function to output the data
def data_output(output_fn, super_dict):
   output = open(output_fn, 'a')
   for codon in super_dict:
      codon_dict = super_dict[codon]
      rev_dict = super_dict[super_dict[codon]["rev_codon"]]
      if abs(float(codon_dict["2Ns"]) + float(rev_dict["2Ns"])) <=0.3:
         closeness = True
      else:
         closeness = False
      data_list = [codon, str(codon_dict["2Ns"]), str(codon_dict["ratio"]), str(codon_dict["#_bin"]), str(codon_dict["total"]), super_dict[codon]["rev_codon"], str(rev_dict["2Ns"]), str(rev_dict["ratio"]), str(rev_dict["#_bin"]), str(rev_dict["total"]), " ", codon,str(codon_dict["2Ns"]), super_dict[codon]["rev_codon"], str(rev_dict["2Ns"]), str(closeness)]
      line = "\t".join(data_list) + "\n"
      #  print(line)
      output.write(line)
   
# Function to analyze the results
def results_analyze(result_file_names, result_dir_path, sfs_dir_path, output_fn):
  output_prep(output_fn)
  super_dict = {}
  for input_file in result_file_names:
    print(input_file)
    input_file_path = result_dir_path + "/" + input_file
    file = open(input_file_path,"r")
    codon = input_file.split(sep="_")[1]
    lines = file.readlines()
    expectation = lines[43].strip().split(sep="\t")
    if expectation[0] == "2Ns":
       pass
    else:
       print(lines[43])
       break
    expectation = expectation[1]
    if abs(float(expectation)) <= 5.0:
       exp_bool = True
    else:
       exp_bool = False
    if lines[49].strip().split(sep="\t")[0] != "i":
       print(lines[49])
       break
    sum = 0.0
    size = 0
    for index in range(50,65):
       data_line = lines[index].strip().split("\t")
       data_ratio = data_line[3]
       fit_ratio = data_line[6]
       if data_ratio == "inf" or fit_ratio == "inf":
          continue 
       data_ratio = float(data_ratio)
       fit_ratio = float(fit_ratio)
       if data_ratio == 0.0 or fit_ratio == 0.0:
          continue
       ratio = 1.0 + (data_ratio - fit_ratio)/data_ratio
       sum += ratio
       size +=1
       pass
    ratio = sum/size
    if ratio <= 1.15 and ratio >= 0.85:
       ratio_bool = True
    else:
       ratio_bool = False
    # data_output(output_fn, codon, expectation, exp_bool, ratio, ratio_bool)
    sfs_file = sfs_dir_path + f"/uk10k_vep_pickorder_flanking_syn_codon_{codon}_down_15"
    sfs_file = open(sfs_file, 'r')
    lines = sfs_file.readlines()
    int_sfs = lines[1].strip().split(sep=" ")
    i = 0
    bin_num = "N/A"
    for bin in int_sfs:
        bin = float(bin)
        if i == 0:
          i+=1
          continue
        elif bin<10:
          bin_num = i
          i = 0
          break
        else:
          i+=1
          pass
    sum = 0
    for bin in int_sfs:
         bin = float(bin)
         if i == 0:
            i+=1
            continue
         else:
            sum += bin
            pass
    rev_codon = codon.split(sep="-")[1] + "-" + codon.split(sep="-")[0]
    super_dict[codon] = {}
    super_dict[codon]["2Ns"] = expectation
    super_dict[codon]["ratio"] = ratio
    super_dict[codon]["#_bin"] = bin_num
    super_dict[codon]["total"] = sum
    super_dict[codon]["rev_codon"] = rev_codon
    pass
  data_output(output_fn, super_dict)

# Function to generate SFS
def generate_sfs(data_fn, output_fn, nc, bank):
   data = open(data_fn, 'r')
   codon_list = []
   syn_sfs = np.zeros(nc, dtype=np.float32)
   int_sfs = np.zeros(nc, dtype=np.float32)
   while True:
      line = data.readline().strip().split(sep="\t")
      if line[0] == "":
         break
      if line[0] == "codon":
         continue
      elif str(line[15]) == "True":
         codon = line[0]
         if codon in codon_list:
            continue
         sfs_fn = f"Document/SFS/FlankingBases/syn_codon/single/{bank}/{bank}_vep_pickorder_flanking_syn_codon_{codon}_down_15"
         sfs = open(sfs_fn, 'r')
         lines = sfs.readlines()
         inter = lines[1].strip().split(sep=" ")
         inter = [float(x) for x in inter]
         syn = lines[3].strip().split(sep=" ")
         syn = [float(x) for x in syn]
         int_sfs += inter
         syn_sfs += syn
        #  codon_list.append(codon)
        #  rev_codon = codon.split(sep="-")
        #  rev_codon = f"{rev_codon[1]}-{rev_codon[0]}"
        #  sfs_fn = f"Document/SFS/FlankingBases/syn_codon/single/{bank}/{bank}_vep_pickorder_flanking_syn_codon_{rev_codon}_down_15"
        #  sfs = open(sfs_fn, 'r')
        #  lines = sfs.readlines()
        #  inter = lines[1].strip().split(sep=" ")
        #  inter = [float(x) for x in inter]
        #  syn = lines[3].strip().split(sep=" ")
        #  syn = [float(x) for x in syn]
        #  int_sfs += inter
        #  syn_sfs += syn
        #  codon_list.append(rev_codon)
         pass
      else:
         continue
   data.close()
   sfs.close()
   output = open(output_fn, 'w')
   int_sfs = [str(x) for x in int_sfs]
   int_sfs = " ".join(int_sfs) + "\n"
   syn_sfs = [str(x) for x in syn_sfs]
   syn_sfs = " ".join(syn_sfs) + "\n"
   output.write("intergenic_variant\n")
   output.write(int_sfs)
   output.write("\nsynonymous_variant\n")
   output.write(syn_sfs)
   output.close()

bank = "uk10k"

   
# result_dir_path = 'Document/SF_Ratios/FlankingBases/syn_codon/uk10k'  # Replace with your directory path
# sfs_dir_path = 'Document/SFS/FlankingBases/syn_codon/single/uk10k'
# result_file_names = list_documents(result_dir_path)
# output_fn = "uk10k_SF_Ratios_data.tsv"
# # print(file_names)
# results_analyze(result_file_names, result_dir_path, sfs_dir_path, output_fn)

# main function
data_fn = f"Document/SF_Ratios/FlankingBases/syn_codon/{bank}_SF_Ratios_data.tsv"
output_sfs = f"{bank}_SF_Ratio_comparison_SFS.txt"
nc = 15
generate_sfs(data_fn, output_sfs, nc, bank)