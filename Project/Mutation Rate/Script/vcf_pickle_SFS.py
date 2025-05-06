# this script porduces the SFS file from pickle file that contains specific factors and its values.
# you can change some assignments like bank, nc,...
# You can also generate SFS with downsampling.
# Mainly two options. Generate with GC type or Codon Pair as primary category

import pickle
import numpy as np
import os 
import logging
from scipy.stats import hypergeom

# assign all the documents in one path, usfel for going through files of each chromosome or other cases. 
def list_documents(directory, extension):
    try:
        file_names = os.listdir(directory)
        file_names = sorted(file_names)
        document_extensions = {extension}
        documents = [f for f in file_names if os.path.isfile(os.path.join(directory, f)) and os.path.splitext(f)[1] in document_extensions]
        return documents
    except Exception as e:
        logging.error(f"Error listing documents: {e}", exc_info=True)

# saving two sfs in text file. SFS in parameter have to be a list of number in each bin.
def sfs_save(output_fn, neutral_sfs, selected_sfs):
   neutral_sfs = [str(x) for x in neutral_sfs]
   neutral_sfs = " ".join(neutral_sfs) + "\n"
   selected_sfs = [str(x) for x in selected_sfs]
   selected_sfs = " ".join(selected_sfs) + "\n"
   with open(output_fn, "w") as f:      
      f.write("synonymous variant\n")
      f.write(selected_sfs)
      f.write("intergenic variant\n")
      f.write(neutral_sfs)

# If you make sfs for each category like GC Type A, B, and n, and you want to combine them.
def sfs_combine_save(directory, output_fn, extension):
   documents = list_documents(directory, extension)
   output = open(output_fn, 'w')
   for input_fn in documents:
      input_fn = os.path.join(directory, input_fn)
      with open(input_fn, "r") as f:
         lines = f.readlines()
         codon_pair = input_fn.split(sep="_")[3].split(sep=".")[0]
      output.write(codon_pair + "\n")
      for line in lines:
        output.write(line)
      output.write("\n")
   output.close()

# Calculate average distance between positions of selected and neutral variants and print it. 
def average_distance_print(distance, count):
   average_distance = distance/count
   print(f"total:{distance}, counts:{count}, average distance: {average_distance}")

# downsampling function
def subsample(pn, obs, nc):
  """
    pn is AN # ,  e.g. 5008
    obs is the AC #,
    nc is the subsample size
  """
  if not (0 < obs < pn):
      return None
  if not (0 < nc <= pn):
      raise ValueError("nc must satisfy 0 < nc <= pn")
  return np.random.hypergeometric(obs, pn - obs, nc)


# One way of pairing up selected variant and neutral variant one by one
def find_and_remove_closest(target, numbers, distance):
    """
    Finds the integer in the list that is closest to the target integer and removes it from the list.
    
    :param target: The target integer
    :param numbers: A list of integers
    :return: The closest integer to the target and the modified list
    """
    if not numbers:
        raise ValueError("The list of numbers is empty")
    
    closest = min(numbers, key=lambda x: abs(x - target))
    numbers.remove(closest)

    distance += abs(target-closest)

    return closest, numbers, distance


# One way of pairing selected and neutral variants with a theoretically minimum distance. 
# Generally better method than the function find_and_remove_closest()
def find_unique_closest_pairs(pos1, pos2):
    """
    Find unique pairs of closest positions between two sorted lists of positions
    using NumPy for efficient distance calculations.
    Returns list of tuples (value_from_pos1, value_from_pos2, distance)
    """

    # Convert lists to NumPy arrays
    pos1 = np.array(pos1)
    pos2 = np.array(pos2)
   
    # Build distance matrix using broadcasting
    distances = np.abs(pos1[:, np.newaxis] - pos2)
   
    pairs = []
    mask = np.ones_like(distances, dtype=bool)
   
    converted_pairs = []

    while True:
        # Find minimum distance in remaining valid positions
        if not np.any(mask):
            break
           
        # Get indices of minimum value in masked array
        min_idx = np.argmin(np.where(mask, distances, np.inf))
        i, j = np.unravel_index(min_idx, distances.shape)
       
        # If no valid pairs left, break
        if distances[i, j] == np.inf:
            break
           
        # Add original values and distance instead of indices
        pairs.append((pos1[i], pos2[j], distances[i, j]))
       
        # Update mask
        mask[i, :] = False  # mask values in row to False
        mask[:, j] = False  # mask values in column to False

        converted_pairs = [(int(a), int(b), int(c)) for a, b, c in pairs]

    return converted_pairs

# Because the find_and_remove_closest function returns a tuple of pairs, this function will give you the position of paired neutral variant corresponding to a selected variant. 
def find_tuple_with_integer(target, tuples_list):
    """
    Finds and returns the first tuple that contains the target integer.
    
    :param target: The integer to find in the list of tuples
    :param tuples_list: The list of tuples to search through
    :return: The first tuple containing the target integer, or None if not found
    """
    for tup in tuples_list:
        if target in tup:
           closest = tup[1]
           distance = tup[2]
           return closest, distance
    return None

# generate list of all synonymous pair
# codon_pair = NNR/NNA or RNN/ANN
# it requires an input text file which has synonymous codon pair in each line. 
def codon_pair_list_generate():
    pair_file_name = "Document/SFS/SynonymousCodon/synonymous_codon_pairs.txt"
    pair_file = open(pair_file_name, 'r')
    lines = pair_file.readlines()
    codon_pair_list = []
    for line in lines:
        line = line.strip().split(sep="\t")
        codon_pair_1 = f'{line[0]}/{line[1]}'
        codon_pair_2 = f'{line[1]}/{line[0]}'
        codon_pair_list.append(codon_pair_1)
        codon_pair_list.append(codon_pair_2)
    return codon_pair_list

# assign CG Type for each synonymous codon pair
def codon_pair_GC_type_dict():
    pair_file_name = "Document/SFS/SynonymousCodon/synonymous_codon_pairs.txt"
    pair_file = open(pair_file_name, 'r')
    lines = pair_file.readlines()
    codon_pair_dict = {}
    for line in lines:
        line = line.strip().split(sep="\t")
        codon_pair_1 = f'{line[0]}/{line[1]}'
        codon_pair_2 = f'{line[1]}/{line[0]}'
        for i in range(3):
          if line[0][i] != line[1][i]:
              GC_type_1 = GC_type_test(line[0][i], line[1][i])
              GC_type_2 = GC_type_test(line[1][i], line[0][i])
              break
        codon_pair_dict[codon_pair_1] = GC_type_1
        codon_pair_dict[codon_pair_2] = GC_type_2
    return codon_pair_dict

# based on ref and alt from vcf file, assign GC Type. A is GC1, B is GC2, n is GC3
def GC_type_test(ref, alt):
  gc = ["G", "C"]
  at = ["A", "T"]
  if ref in gc and alt in gc or ref in at and alt in at:
      return "n"
  elif ref in at and alt in gc:
      return "A"
  elif ref in gc and alt in at:
      return "B"

# By using an interval dictionary that can be obtained from mask file, you can eliminate the variants whose positions are not in mask file.
def bed_interval_eliminate(interval_dict, pos_list):
   for pos in pos_list:
      try:
        pos_bool = interval_dict[pos]
        if pos_bool:
           continue
        else:
           pos_list.remove(pos)
           pass
      except KeyError:
        continue
   return pos_list

# Generate SFS based on GC bias type, A, B or n
#   Pick synonymous_variant
#     Pick intergenic_variant with same GC bias type, mutation rate bin number, and closest pos
#       Generate SFS for syn and int
def SFS_GC_bias_type(syn_dict, int_dict, syn_sfs, int_sfs, GC_bias, pairing_method, distance, count, interval_dict):
  try:
      for codon_pair in syn_dict[GC_bias]:
         for syn_mut in syn_dict[GC_bias][codon_pair]:
          int_pos_list = list(int_dict[GC_bias][syn_mut].keys())
          int_pos_list = bed_interval_eliminate(interval_dict, int_pos_list)
          if pairing_method == "one by one":
            for syn_pos in syn_dict[GC_bias][codon_pair][syn_mut]:
                syn_an_ac = syn_dict[GC_bias][codon_pair][syn_mut][syn_pos].split(sep="_")
                syn_ac = int(syn_an_ac[1])
                int_pos, int_pos_list, distance = find_and_remove_closest(syn_pos, int_pos_list, distance)
                int_an_ac = int_dict[GC_bias][syn_mut][int_pos].split(sep="_")
                int_ac = int(int_an_ac[1])
                syn_sfs[syn_ac] +=1
                int_sfs[int_ac] +=1
                count += 1
                pass
          elif pairing_method == "average":
             syn_pos_list = list(syn_dict[GC_bias][codon_pair][syn_mut].keys())
             syn_pos_list = bed_interval_eliminate(interval_dict, syn_pos_list)
            #  print(syn_pos_list[0:20])
            #  print(int_pos_list[0:20])
             pairs = find_unique_closest_pairs(syn_pos_list, int_pos_list)
             for syn_pos in syn_pos_list:
                int_pos, individual_distance = find_tuple_with_integer(syn_pos, pairs)
                syn_an_ac = syn_dict[GC_bias][codon_pair][syn_mut][syn_pos].split(sep="_")
                syn_ac = int(syn_an_ac[1])
                int_an_ac = int_dict[GC_bias][syn_mut][int_pos].split(sep="_")
                int_ac = int(int_an_ac[1])
                syn_sfs[syn_ac] +=1
                int_sfs[int_ac] +=1
                count += 1
                distance += individual_distance
                pass
          else:
             print(f"Unknown pairing method : {pairing_method}")
      return syn_sfs, int_sfs, distance, count
  except KeyError:
      return syn_sfs, int_sfs, distance, count
  


#   Pick synonymous_variant
#     Pick intergenic_variant with same GC bias type, mutation rate bin number, and closest pos
#       Generate SFS for syn and int with downsampling
def SFS_GC_bias_type_downsample(syn_dict, int_dict, syn_sfs, int_sfs, GC_bias, pairing_method, distance, count, interval_dict, nc):
  try:
      for codon_pair in syn_dict[GC_bias]:
         for syn_mut in syn_dict[GC_bias][codon_pair]:
          int_pos_list = list(int_dict[GC_bias][syn_mut].keys())
          int_pos_list = bed_interval_eliminate(interval_dict, int_pos_list)
          if pairing_method == "one by one":
            for syn_pos in syn_dict[GC_bias][codon_pair][syn_mut]:
                syn_an_ac = syn_dict[GC_bias][codon_pair][syn_mut][syn_pos].split(sep="_")
                syn_ac = int(syn_an_ac[1])
                int_pos, int_pos_list, distance = find_and_remove_closest(syn_pos, int_pos_list, distance)
                int_an_ac = int_dict[GC_bias][syn_mut][int_pos].split(sep="_")
                int_ac = int(int_an_ac[1])
                an = int(syn_an_ac[0])
                syn_ac = subsample(an, syn_ac, nc)
                int_ac = subsample(an, int_ac, nc)
                if syn_ac == None or int_ac == None:
                   continue
                syn_sfs[syn_ac] +=1
                int_sfs[int_ac] +=1
                count += 1
                pass
          elif pairing_method == "average":
             syn_pos_list = list(syn_dict[GC_bias][codon_pair][syn_mut].keys())
             syn_pos_list = bed_interval_eliminate(interval_dict, syn_pos_list)
            #  print(syn_pos_list[0:20])
            #  print(int_pos_list[0:20])
             pairs = find_unique_closest_pairs(syn_pos_list, int_pos_list)
             for syn_pos in syn_pos_list:
                int_pos, individual_distance = find_tuple_with_integer(syn_pos, pairs)
                syn_an_ac = syn_dict[GC_bias][codon_pair][syn_mut][syn_pos].split(sep="_")
                syn_ac = int(syn_an_ac[1])
                int_an_ac = int_dict[GC_bias][syn_mut][int_pos].split(sep="_")
                int_ac = int(int_an_ac[1])
                an = int(syn_an_ac[0])
                syn_ac = subsample(an, syn_ac, nc)
                int_ac = subsample(an, int_ac, nc)
                if syn_ac == None or int_ac == None:
                   continue
                syn_sfs[syn_ac] +=1
                int_sfs[int_ac] +=1
                count += 1
                distance += individual_distance
                pass
          else:
             print(f"Unknown pairing method : {pairing_method}")
      return syn_sfs, int_sfs, distance, count
  except KeyError:
      return syn_sfs, int_sfs, distance, count


# Generate SFS based on codon_pair
#   For each codon pair
#     Pick synonymous_variant
#       Pick intergenic_variant with same GC bias type, mutation rate bin number, and closest pos
#         Generate SFS for syn and int
def SFS_codon_pair_GC_bias(syn_dict, int_dict, syn_sfs, int_sfs, codon_pair, GC_bias, pairing_method, distance, count, interval_dict):
    try:
      for syn_mut in syn_dict[GC_bias][codon_pair]:
         int_pos_list = list(int_dict[GC_bias][syn_mut].keys())
         int_pos_list = bed_interval_eliminate(interval_dict, int_pos_list)
         if pairing_method == "one by one":
           for syn_pos in syn_dict[GC_bias][codon_pair][syn_mut]:
              syn_an_ac = syn_dict[GC_bias][codon_pair][syn_mut][syn_pos].split(sep="_")
              syn_ac = int(syn_an_ac[1])
              int_pos, int_pos_list, distance = find_and_remove_closest(syn_pos, int_pos_list, distance)
              int_an_ac = int_dict[GC_bias][syn_mut][int_pos].split(sep="_")
              int_ac = int(int_an_ac[1])
              syn_sfs[syn_ac] +=1
              int_sfs[int_ac] +=1
              count += 1
              pass
         elif pairing_method == "average":
            syn_pos_list = list(syn_dict[GC_bias][codon_pair][syn_mut].keys())
            syn_pos_list = bed_interval_eliminate(interval_dict, syn_pos_list)
            # print(syn_pos_list[0:20])
            # print(int_pos_list[0:20])
            pairs = find_unique_closest_pairs(syn_pos_list, int_pos_list)
            for syn_pos in syn_pos_list:
                int_pos, individual_distance = find_tuple_with_integer(syn_pos, pairs)
                syn_an_ac = syn_dict[GC_bias][codon_pair][syn_mut][syn_pos].split(sep="_")
                syn_ac = int(syn_an_ac[1])
                int_an_ac = int_dict[GC_bias][syn_mut][int_pos].split(sep="_")
                int_ac = int(int_an_ac[1])
                syn_sfs[syn_ac] +=1
                int_sfs[int_ac] +=1
                count += 1
                distance += individual_distance
                pass
         else: 
            print(f"Unknown pairing method: {pairing_method}")
      return syn_sfs, int_sfs, distance, count
    except KeyError:
      return syn_sfs, int_sfs, distance, count
   
# Generate SFS based on codon_pair
#   For each codon pair
#     Pick synonymous_variant
#       Pick intergenic_variant with same GC bias type, mutation rate bin number, and closest pos
#         Generate SFS for syn and int with downsampling
def SFS_codon_pair_GC_bias_downsample(syn_dict, int_dict, syn_sfs, int_sfs, codon_pair, GC_bias, pairing_method, distance, count, interval_dict, nc):
    try:
      for syn_mut in syn_dict[GC_bias][codon_pair]:
         int_pos_list = list(int_dict[GC_bias][syn_mut].keys())
         int_pos_list = bed_interval_eliminate(interval_dict, int_pos_list)
         if pairing_method == "one by one":
           for syn_pos in syn_dict[GC_bias][codon_pair][syn_mut]:
              syn_an_ac = syn_dict[GC_bias][codon_pair][syn_mut][syn_pos].split(sep="_")
              syn_ac = int(syn_an_ac[1])
              int_pos, int_pos_list, distance = find_and_remove_closest(syn_pos, int_pos_list, distance)
              int_an_ac = int_dict[GC_bias][syn_mut][int_pos].split(sep="_")
              int_ac = int(int_an_ac[1])
              an = int(syn_an_ac[0])
              syn_ac = subsample(an, syn_ac, nc)
              int_ac = subsample(an, int_ac, nc)
              if syn_ac == None or int_ac == None:
                   continue
              syn_sfs[syn_ac] +=1
              int_sfs[int_ac] +=1
              count += 1
              pass
         elif pairing_method == "average":
            syn_pos_list = list(syn_dict[GC_bias][codon_pair][syn_mut].keys())
            syn_pos_list = bed_interval_eliminate(interval_dict, syn_pos_list)
            # print(syn_pos_list[0:20])
            # print(int_pos_list[0:20])
            pairs = find_unique_closest_pairs(syn_pos_list, int_pos_list)
            for syn_pos in syn_pos_list:
                int_pos, individual_distance = find_tuple_with_integer(syn_pos, pairs)
                syn_an_ac = syn_dict[GC_bias][codon_pair][syn_mut][syn_pos].split(sep="_")
                syn_ac = int(syn_an_ac[1])
                int_an_ac = int_dict[GC_bias][syn_mut][int_pos].split(sep="_")
                int_ac = int(int_an_ac[1])
                an = int(syn_an_ac[0])
                syn_ac = subsample(an, syn_ac, nc)
                int_ac = subsample(an, int_ac, nc)
                if syn_ac == None or int_ac == None:
                   continue
                syn_sfs[syn_ac] +=1
                int_sfs[int_ac] +=1
                count += 1
                distance += individual_distance
                pass
         else: 
            print(f"Unknown pairing method: {pairing_method}")
      return syn_sfs, int_sfs, distance, count
    except KeyError:
      return syn_sfs, int_sfs, distance, count

# If you want to generate SFS based on GC Type, run this function
def SFS_GC_type_script():
  bank = "1kG"
  nc = 200
  pickle_directory = "Document/Pickle/MutationRate/SuperDictionary/1kG_GC_syn_int"
  documents = list_documents(pickle_directory, ".p")
  interval_pickle = "Document/Pickle/MutationRate/PosDict/1kG_super_dict_bed_pilot.p"
  # assign pairing method to "one by one" or "average"
  pairing_method = "average"
  GC_type_list = ["A", "B", "n"]
  for GC_type in GC_type_list:
    print(f"GC_type:{GC_type}")
    syn_sfs = [0]*(nc+1)
    int_sfs = [0]*(nc+1)
    distance = 0
    count = 0
    for chr in range(1,23):
        std_chr = f"chr{chr}"
        with open(interval_pickle, 'rb') as f:
         interval_dict = pickle.load(f)[std_chr]
        for input_fn in documents:
            # print(input_fn)
            if input_fn.split(sep="_")[1] != std_chr:
              continue
            elif input_fn.split(sep="_")[3] == "syn.p":
              input_fn = os.path.join(pickle_directory, input_fn)
              with open(input_fn, "rb") as f:
                  syn_dict = pickle.load(f)
            elif input_fn.split(sep="_")[3] == "int.p":
              input_fn = os.path.join(pickle_directory, input_fn)
              with open(input_fn, "rb") as f:
                  int_dict = pickle.load(f)
            else:
              print(f"ERROR: unknown file type, {input_fn}")
        # syn_sfs, int_sfs, distance, count = SFS_GC_bias_type(syn_dict, int_dict, syn_sfs, int_sfs, GC_type, pairing_method, distance, count, interval_dict)
        syn_sfs, int_sfs, distance, count = SFS_GC_bias_type_downsample(syn_dict, int_dict, syn_sfs, int_sfs, GC_type, pairing_method, distance, count, interval_dict, nc)
    # average_distance_print(distance, count)
    output_fn = f"Document/SFS/MutationRate/GCBias/Mask/Pilot/{bank}_GC_bias_downsample_200_{GC_type}.txt"
    sfs_save(output_fn, int_sfs, syn_sfs)
    print("end")

  directory = "Document/SFS/MutationRate/GCBias/Mask/Pilot"
  output_file = "1kG_GC_bias_type_pilot_downsample_200.txt"
  sfs_combine_save(directory, output_file, ".txt")

# If you want to generate SFS based on synonymous codon pair, run this function.
def SFS_codon_pair_GC_script():
  bank = "1kG"
  nc = 200
  pickle_directory = "Document/Pickle/MutationRate/SuperDictionary/1kG_GC_syn_int"
  documents = list_documents(pickle_directory, ".p")
  # assign pairing method to "one by one" or "average"
  interval_pickle = "Document/Pickle/MutationRate/PosDict/1kG_super_dict_bed_pilot.p"
  pairing_method = "average"
  codon_pair_list = codon_pair_list_generate()
  codon_pair_dict = codon_pair_GC_type_dict()
  for codon_pair in codon_pair_list:
    print(f"{codon_pair} start")
    syn_sfs = [0]*(nc+1)
    int_sfs = [0]*(nc+1)
    distance = 0
    count = 0
    GC_bias = codon_pair_dict[codon_pair]
    for chr in range(1, 23):
      std_chr = f"chr{chr}"
      with open(interval_pickle, 'rb') as f:
         interval_dict = pickle.load(f)[std_chr]
      for input_fn in documents:
          # print(input_fn)
          if input_fn.split(sep="_")[1] != std_chr:
            continue
          elif input_fn.split(sep="_")[3] == "syn.p":
            input_fn = os.path.join(pickle_directory, input_fn)
            with open(input_fn, "rb") as f:
                syn_dict = pickle.load(f)
          elif input_fn.split(sep="_")[3] == "int.p":
            input_fn = os.path.join(pickle_directory, input_fn)
            with open(input_fn, "rb") as f:
                int_dict = pickle.load(f)
          else:
            print(f"ERROR: unknown file type, {input_fn}")
      # syn_sfs, int_sfs, distance, count = SFS_codon_pair_GC_bias(syn_dict, int_dict, syn_sfs, int_sfs, codon_pair, GC_bias, pairing_method, distance, count, interval_dict)
      syn_sfs, int_sfs, distance, count = SFS_codon_pair_GC_bias_downsample(syn_dict, int_dict, syn_sfs, int_sfs, codon_pair, GC_bias, pairing_method, distance, count, interval_dict, nc)
    # average_distance_print(distance, count)
    codon_pair = codon_pair.replace("/", "-")
    output_fn = f"Document/SFS/MutationRate/CodonPair/Mask/Pilot/{bank}_GC_bias_downsample_200_{codon_pair}.txt"
    sfs_save(output_fn, int_sfs, syn_sfs)
    print("end")
  
  directory = "Document/SFS/MutationRate/CodonPair/Mask/Pilot"
  output_file = "1kG_GC_codon_pair_pilot_downsample_200.txt"
  sfs_combine_save(directory, output_file, ".txt")



# run one of below, or both
# SFS_GC_type_script()
# SFS_codon_pair_GC_script()