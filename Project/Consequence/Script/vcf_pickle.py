# this script reads the binary gzipped file, creats a super dictionary where the first keys are consequences, the second keys are AC and AN, and values are the frequency, and saves the dictionary to pickle file.
# The script also creates a super dictionary where the first keys are GC bias, the second keys are codon pairs, the third keys are mutation rate bin numbers, and the values are positions and AC/AN.

import gzip
import pysam
import pickle
import os
import random
import paramiko
import io
import logging
import gc
 

# make a pickle file containing a super dictionary with CONSEQUENCE as first key and AC_AN as second key.
def conseq_pickle(output_fn):
    zero_one = [0, 1]
    super_dict = {}
    skipped = 0
    for i in range(1, 23):
        print("start ",i,end="")
        fn1 = os.path.join("vcf/","1kG-chr{}-vep-every.vcf.gz".format(i))
        gunzipped = gzip.open(fn1, 'r')
        f1 = pysam.VariantFile(gunzipped)
        while True:
            try:
                record = next(f1)
            except OSError:
                continue
            except StopIteration:
                break
            try:
                vep = record.info["CSQ"]
                # vep = record.info["vep"]
            except KeyError:
                continue
            # if len(vep) > 1: # skip all records that indicate more than one function
            #     skipped +=1
            #     continue
            if len(record.ref) != 1 or any(len(alt) != 1 for alt in record.alts):
                continue
            conseq_list = []
            for annotation in vep:
                conseq = annotation.split(sep="|")[1]
                if conseq in conseq_list:
                    continue
                else:
                    conseq_list.append(conseq)
                    ac = record.info["AC"][0]
                    an = record.info["AN"]
                    data = "AC_" + str(ac) + "_AN_" + str(an)
                    if conseq in super_dict:
                        if data in super_dict[conseq]:
                            super_dict[conseq][data] += 1
                            pass
                        else:
                            super_dict[conseq][data] = 1
                            # print(random.choice(zero_one))
                            pass
                    elif conseq not in super_dict:
                        super_dict[conseq] = {}
                        super_dict[conseq][data] = 1
                        pass
        print(f" done, {skipped}")
    pickle.dump(super_dict, open(output_fn, 'wb'))

# make a pickle file containing a super dictionary with only specific CONSEQUENCE as first key and AC_AN as second key.
def conseq_pickle_specific_conseq(output_fn, spe_conseq_list):
    zero_one = [0, 1]
    key_error = 0
    conseq_unmatched = 0
    ac_error = 0
    vep_error = 0
    super_dict = {}
    for spe_conseq in spe_conseq_list:
        super_dict[spe_conseq] = {}
    for i in range(1, 23):
        print("start ",i,end="")
        fn1 = os.path.join("/mnt/shared_dir/vcf/uk10k/uk10k_pantro6filtered_regfiltered/uk10k_chr{}_official_merged_vep_pickorder_fixed_pantro6filtered_regfiltered_modified.vcf.bgz".format(i))
        gunzipped = gzip.open(fn1, 'r')
        f1 = pysam.VariantFile(gunzipped)
        for record in f1:
            try:
                vep = record.info["CSQ"]
                # vep = record.info["vep"]
            except KeyError:
                key_error +=1
                continue
            if len(record.info["AC"])>1:
                ac_error +=1
                continue
            else:
                ac = record.info["AC"][0]
            an = record.info["AN"]
            data = "AC_" + str(ac) + "_AN_" + str(an)
            if len(vep)>1:
                vep_error +=1
                continue
            else:
                conseq = vep[0].split(sep="|")[1]
            if conseq not in spe_conseq_list:
                conseq_unmatched +=1
                continue
            else:
                if data in super_dict[conseq]:
                    super_dict[conseq][data] += 1
                    pass
                else:
                    super_dict[conseq][data] = 1
                    print(random.choice(zero_one))
                    pass
                   
        print(f" done, key_error: {key_error}, ac_error: {ac_error}, vep_error: {vep_error}, unmatched: {conseq_unmatched}")
    pickle.dump(super_dict, open(output_fn, 'wb'))

# make a pickle file containing a super dictionary with CONSEQUENCE as first key and FOLDED AC_AN as second key.
def conseq_pickle_folded(output_fn):
    super_dict = {}
    for i in range(1, 23):
        print("start ",i,end="")
        fn1 = os.path.join(f"vcf/uk10k_chr{i}_official_merged_vep_pickorder_fixed_pantro6filtered_regfiltered_modified.vcf.bgz")
        gunzipped = gzip.open(fn1, 'r')
        f1 = pysam.VariantFile(gunzipped)
        for record in f1:
            try:
                vep = record.info["CSQ"]
                # vep = record.info["vep"]
            except KeyError:
                continue
            if len(vep) > 1: # skip all records that indicate more than one function
                continue
            annotation = vep[0]
            conseq = annotation.split(sep="|")[1]
            if len(record.info["AC"]) == 1:
                ac = record.info["AC"][0]
            else: 
                continue
            an = record.info["AN"]/2
            if ac > an:
                ac = record.info["AN"] - record.info["AC"][0]
            # print(mut_seq)
            if conseq in super_dict:
                data = "AC_" + str(ac) + "_AN_" + str(an)
                if data in super_dict[conseq]:
                    super_dict[conseq][data] += 1
                    pass
                else:
                    super_dict[conseq][data] = 1
                    pass
            if conseq not in super_dict:
                data = "AC_" + str(ac) + "_AN_" + str(an)
                super_dict[conseq] = {}
                super_dict[conseq][data] = 1
                pass
        print(" done")
    pickle.dump(super_dict, open(output_fn, 'wb'))

def conseq_pickle_sift4g(output_pickle):
  zero_one = [0, 1]
  super_dict = {}
  for i in range(1, 23):
    print("start ",i,end="")
    fn1 = f"vcf/1kG-chr{i}-sift4g.vcf.gz"
    gunzzipped = gzip.open(fn1, "r")
    vcf = pysam.VariantFile(gunzzipped)
    # print(vcf.header)  
    skipped = 0
    for record in vcf:
      try:
        info = record.info["SIFTINFO"][0]
        # print(info)
      except KeyError:
        skipped +=1
        continue
      conseq = info.split(sep="|")[5]
      ac = record.info["AC"][0]
      an = record.info["AN"]
      # print(mut_seq)
      data = "AC_" + str(ac) + "_AN_" + str(an)
      if conseq in super_dict:
          if data in super_dict[conseq]:
              super_dict[conseq][data] += 1
              pass
          else:
              super_dict[conseq][data] = 1
              print(random.choice(zero_one))
              pass
      if conseq not in super_dict:
          super_dict[conseq] = {}
          super_dict[conseq][data] = 1
          pass
    vcf.close()
    print(" done")
    print(f"skipped {skipped}times")
  pickle.dump(super_dict, open(output_pickle, 'wb'))


out_pickle = "1kG-sift4g-single-conseq.p"
conseq_pickle_sift4g(out_pickle)      