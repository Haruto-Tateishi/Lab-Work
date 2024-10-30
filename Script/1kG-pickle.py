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
    skipped = 0
    for i in range(1, 23):
        print("start ",i,end="")
        fn1 = op.join("vcf/","1kG-chr{}-vep-every.vcf.gz".format(i))
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
                continueac = record.info["AC"][0]
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
                if conseq not in super_dict:
                    super_dict[conseq] = {}
                    super_dict[conseq][data] = 1
                   
            conseq_list = []        # While loop
            if len(record.info["AC"]) == 1:
                ac = record.info["AC"][0]
            else: 
                continue
            for annotation in vep:
                    conseq = annotation.split(sep="|")[1]
                    conseq_list.append(conseq)
            for spe_conseq in spe_conseq_list:
                if spe_conseq not in conseq_list:
                    continue
                else:
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

def conseq_flanking_pickle(output_fn):
    zero_one = [0, 1]
    super_dict = {}
    for chr in range(1, 23):
        print(f"start{chr} \n")
        position_skipped =0
        format_skipped=0
        key_error_skipped=0
        sequence_skipped = 0
        fn1 = op.join("vcf/1kG_chr{}_vep_pickorder_pantro6filtered_regfiltered_modified_codon.vcf.gz".format(chr))
        gunzipped = gzip.open(fn1, 'r')
        f1 = pysam.VariantFile(gunzipped)
        super_dict[chr] = {}
        flanking_list = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]
        variant_list = ["A/C", "A/G", "A/T", "C/A", "C/G", "C/T", "G/A", "G/C", "G/T", "T/A", "T/C", "T/G"]
        for variant in variant_list:
            for bases in flanking_list:
                flanking = f"{bases[0]}{variant}{bases[1]}"
                super_dict[chr][flanking] = {}
        skip_list = []
        pre_pos = 0
        while True:
            try:
                record = next(f1)
            except OSError:
                continue
            except StopIteration:
                break
            cur_pos = record.pos
            if cur_pos - pre_pos == 1:
                skip_list.append(pre_pos)
                skip_list.append(cur_pos)
                # print("prereading")
            pre_pos = cur_pos
            pass
        print("prework done")
        f1.close()
        gunzipped = None
        gunzipped = gzip.open(fn1, 'r')
        f1 = pysam.VariantFile(gunzipped)
        while True:
            try:
                record = next(f1)
            except OSError:
                continue
            except StopIteration:
                break
            if record.pos in skip_list:
                position_skipped+=1
                continue
            try:
                vep = record.info["CSQ"]
                # vep = record.info["vep"]
            except KeyError:
                continue
            if len(vep) > 1: # skip all records that indicate more than one function
                format_skipped +=1
                continue
            if len(record.ref) != 1 or any(len(alt) != 1 for alt in record.alts):
                format_skipped +=1
                continue
            annotation = vep[0]
            try:
                conseq = annotation.split(sep="|")[1]
                ac = record.info["AC"][0]
                pos = record.pos
                refseq = record.info["RefSeq"].upper()
                # print(refseq)
                if "N" in refseq:
                    sequence_skipped+=1
                    break
                ref = record.ref
                alt = record.alts[0]
            except KeyError:
                    key_error_skipped+=1
                    continue
            flanking = f"{refseq[0]}{ref}/{alt}{refseq[2]}"
            # print(f"ref:{ref} alt: {alt}")
            if conseq in super_dict[chr][flanking]:
                super_dict[chr][flanking][conseq][pos] = ac
                continue
            elif conseq not in super_dict[chr][flanking]:
                super_dict[chr][flanking][conseq] = {}
                super_dict[chr][flanking][conseq][pos] = ac
                continue
        print(f" done, skipped, position:{position_skipped}, format:{format_skipped}, sequence:{sequence_skipped}, key_error:{key_error_skipped}")
    pickle.dump(super_dict, open(output_fn, 'wb'))


output = "1kG_vep_pickorder_pantro6filtered_flanking_location.p"
conseq_flanking_pickle(output)