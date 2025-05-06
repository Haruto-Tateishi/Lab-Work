# this script reads the binary gzipped file, creats a super dictionary where the first keys are consequences, the second keys are AC and AN, and values are the frequency, and saves the dictionary to pickle file.
# The script also creates a super dictionary where the first keys are GC bias, the second keys are codon pairs, the third keys are mutation rate bin numbers, and the values are positions and AC/AN.

import gzip
import pysam
import pickle
import os
import random

# make a pickle file containing a super dictionary with SYNONYMOUS CODON PAIR as first key and AC_AN as second key.
def codon_pair_pickle(output_fn):
    pair_file_name = "Document/SFS/SynonymousCodon/synonymous_codon_pairs.txt"
    pair_file = open(pair_file_name, 'r')
    lines = pair_file.readlines()
    super_dict = {}
    for line in lines:
        line = line.strip().split(sep="\t")
        codon_pair_1 = f'{line[0]}->{line[1]}'
        codon_pair_2 = f'{line[1]}->{line[0]}'
        super_dict[codon_pair_1] = {}
        super_dict[codon_pair_2] = {}

    # print(super_dict)
    exception_list = []
    for i in range(1, 23):
        print("start ",i,end="")
        fn1 = os.path.join("/home/haruto/Lab-Work/vcf/","1kG-chr{}-pantro6-filtered-syn.vcf".format(i))
        f1 = pysam.VariantFile(fn1)
        for record in f1:
            try:
                vep = record.info["CSQ"]
                # vep = record.info["vep"]
            except KeyError:
                continue
            # if len(vep) > 1: # skip all records that indicate more than one function
            #     continue
            if len(record.info["AC"]) == 1:
                ac = record.info["AC"][0]
            else: 
                continue
            for annotation in vep:
                conseq = annotation.split(sep="|")[1]
                if conseq != "synonymous_variant":
                    continue
                else:
                    an = record.info["AN"]
                    ref_seq = record.info["RdefefSeq"].upper()
                    # print(ref_seq)
                    mut_seq = record.info["MutSeq"]
                    if type(mut_seq) == tuple:
                        mut_seq = mut_seq[0].split(sep=",")
                    else:
                        mut_seq = mut_seq.split(sep=",")
                    # print(mut_seq)
                    if len(mut_seq) > 1:
                        continue
                    else:
                        mut_seq = mut_seq[0].upper()
                        # print(mut_seq)
                    codon_pair = ref_seq + "->" + mut_seq
                    if codon_pair in super_dict:
                        data = "AC_" + str(ac) + "_AN_" + str(an)
                        if data in super_dict[codon_pair]:
                            super_dict[codon_pair][data] += 1
                        else:
                            super_dict[codon_pair][data] = 1
                    else:
                        if codon_pair not in exception_list:
                            exception_list.append(codon_pair)
                    break
        print(" done")
    print(exception_list)
    pickle.dump(super_dict, open(output_fn, 'wb'))    
