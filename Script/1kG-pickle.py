# this script reads the binary gzipped file, creats a super dictionary where the first keys are consequences, the second keys are AC and AN, and values are the frequency, and saves the dictionary to pickle file.

import gzip
import pysam
import pickle
import os.path as op
import random

def codon_pair_pickle(output_fn):
    pair_file_name = "synonymous_codon_pairs.txt"
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
        fn1 = op.join("/home/haruto/Lab-Work/vcf/","1kG-chr{}-pantro6-filtered-syn.vcf".format(i))
        f1 = pysam.VariantFile(fn1)
        for record in f1:
            try:
                vep = record.info["CSQ"]
                # vep = record.info["vep"]
            except KeyError:
                continue
            # if len(vep) > 1: # skip all records that indicate more than one function
            #     continue
            for annotation in vep:
                conseq = annotation.split(sep="|")[1]
                if conseq != "synonymous_variant":
                    continue
                else:
                    ac = record.info["AC"][0]
                    an = record.info["AN"]
                    ref_seq = record.info["RefSeq"].upper()
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

def conseq_pickle(output_fn):
    zero_one = [0, 1]
    super_dict = {}
    for i in range(1, 23):
        print("start ",i,end="")
        fn1 = op.join("/home/haruto/Lab-Work/vcf/","1kG-chr{}-vep-every.vcf.gz".format(i))
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
        print(" done")
    pickle.dump(super_dict, open(output_fn, 'wb'))

def conseq_pickle_specific_conseq(output_fn, spe_conseq_list):
    zero_one = [0, 1]
    super_dict = {}
    for spe_conseq in spe_conseq_list:
        super_dict[spe_conseq] = {}
    for i in range(1, 23):
        print("start ",i,end="")
        fn1 = op.join("/home/haruto/Lab-Work/vcf/","1kG-chr{}-vep-every.vcf.gz".format(i))
        gunzipped = gzip.open(fn1, 'r')
        f1 = pysam.VariantFile(gunzipped)
        for record in f1:
            try:
                vep = record.info["CSQ"]
                # vep = record.info["vep"]
            except KeyError:
                continue
            conseq_list = []
            for annotation in vep:
                    conseq = annotation.split(sep="|")[1]
                    conseq_list.append(conseq)
            for spe_conseq in spe_conseq_list:
                if spe_conseq not in conseq_list:
                    continue
                else:
                    ac = record.info["AC"][0]
                    an = record.info["AN"]
                    # print(mut_seq)
                    data = "AC_" + str(ac) + "_AN_" + str(an)
                    if data in super_dict[spe_conseq]:
                        super_dict[spe_conseq][data] += 1
                        pass
                    else:
                        super_dict[spe_conseq][data] = 1
                        print(random.choice(zero_one))
                        pass
        print(" done")
    pickle.dump(super_dict, open(output_fn, 'wb'))

def conseq_pickle_folded(output_fn):
    super_dict = {}
    for i in range(1, 23):
        print("start ",i,end="")
        fn1 = op.join("/home/haruto/Lab-Work/vcf/","1kG-chr{}-vep-every.vcf.gz".format(i))
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
            ac = record.info["AC"][0]
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

def conseq_seq_pickle(output_fn):
    zero_one = [0, 1]
    super_dict = {}
    for i in range(1, 23):
        print("start ",i,end="")
        fn1 = op.join("/home/haruto/Lab-Work/vcf/","1kG-chr{}-vep-every.vcf.gz".format(i))
        gunzipped = gzip.open(fn1, 'r')
        f1 = pysam.VariantFile(gunzipped)
        for record in f1:
            try:
                seq = record.info["RefSeq"]
                # vep = record.info["vep"]
            except KeyError:
                continue
            try:
                vep = record.info["CSQ"]
                # vep = record.info["vep"]
            except KeyError:
                continue
            if len(vep) > 1: # skip all records that indicate more than one function
                continue
            if len(seq) > 1: # skip all records that indicate more than one alt allele
                continue
            alt = record.ref
            codon_seq = list(seq)
            codon_seq.insert(2, f"/{alt}")
            codon_seq = "".join(codon_seq)
            annotation = vep[0]
            conseq = annotation.split(sep="|")[1]
            ac = record.info["AC"][0]
            an = record.info["AN"]
            data = "AC_" + str(ac) + "_AN_" + str(an)
            # print(mut_seq)
            if conseq in super_dict:
                if codon_seq in super_dict[conseq]:
                    if data in super_dict[conseq][codon_seq]:
                        super_dict[conseq][codon_seq][data] += 1
                        pass
                    else:
                        super_dict[conseq][codon_seq][data] = 1
                        print(random.choice(zero_one))
                else:
                    super_dict[conseq][codon_seq] = {}
                    super_dict[conseq][codon_seq][data] = 1
                    pass
            else:
                super_dict[conseq] = {}
                super_dict[conseq][codon_seq] = {}
                super_dict[conseq][codon_seq][data] = 1
                pass
        print(" done")
    pickle.dump(super_dict, open(output_fn, 'wb'))


output = "1kG_vep_syn_mis_ANAC_counts_all.p"
spe_conseq_list = ["synonymous_variant", "missense_variant"]
conseq_pickle_specific_conseq(output, spe_conseq_list)