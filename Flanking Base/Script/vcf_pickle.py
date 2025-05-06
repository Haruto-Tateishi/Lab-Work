# this script reads the binary gzipped file, creats a super dictionary where the first keys are consequences, the second keys are AC and AN, and values are the frequency, and saves the dictionary to pickle file.
# The script also creates a super dictionary where the first keys are GC bias, the second keys are codon pairs, the third keys are mutation rate bin numbers, and the values are positions and AC/AN.

import gzip
import pysam
import pickle
import os

# make a pickle file containing a super dictionary with FLANKING BASE as fist key, CONSEQUENCE as second key, position as third key and AC as value.
def conseq_flanking_pickle(output_fn):
    zero_one = [0, 1]
    super_dict = {}
    for chr in range(1, 23):
        print(f"start{chr} \n")
        position_skipped =0
        format_skipped=0
        key_error_skipped=0
        sequence_skipped = 0
        fn1 = os.path.join("/mnt/shared_dir/vcf/uk10k/uk10k_pantro6filtered_regfiltered/uk10k_pantro6filtered_regfiltered_codon/uk10k_chr{}_official_merged_vep_pickorder_fixed_pantro6filtered_regfiltered_modified_codon.vcf.gz".format(chr))
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
        pre_conseq = None
        pre_flank = None
        print("prework done")
        pre_pos = 0
        while True:
            try:
                record = next(f1)
            except OSError:
                continue
            except StopIteration:
                break
            # if record.pos in skip_list:
            #     position_skipped+=1
            #     continue
            try:
                vep = record.info["CSQ"]
                # vep = record.info["vep"]
            except KeyError:
                continue
            cur_pos = record.pos
            if cur_pos-pre_pos<=1:
                pre_pos = cur_pos
                position_skipped +=1
                try:
                    del super_dict[chr][pre_flank][pre_conseq][pre_pos]
                    continue
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
                super_dict[chr][flanking][conseq][cur_pos] = ac
                pre_flank = flanking
                pre_conseq = conseq
                pre_pos = cur_pos
                continue
            elif conseq not in super_dict[chr][flanking]:
                super_dict[chr][flanking][conseq] = {}
                super_dict[chr][flanking][conseq][cur_pos] = ac
                pre_flank = flanking
                pre_conseq = conseq
                pre_pos = cur_pos
                # print(random.choice(zero_one))
                continue
        print(f" done, skipped, position:{position_skipped}, format:{format_skipped}, sequence:{sequence_skipped}, key_error:{key_error_skipped}")
    pickle.dump(super_dict, open(output_fn, 'wb'))

# make a pickle file containing a super dictionary with FLANKING BASE as fist key, CONSEQUENCE as second key, position as third key and AC as value.
# if conseq is synonymous variant, the codon pair is also added to the dictionary
def conseq_flanking_syn_pair_pickle(output_fn, con, exp):
    zero_one = [0, 1]
    super_dict = {}
    for chr in range(1, 23):
        print(f"start{chr} \n")
        con_processed = 0
        exp_processed = 0
        position_skipped =0
        format_skipped=0
        key_error_skipped=0
        sequence_skipped = 0
        conseq_skipped = 0
        fn1 = os.path.join("vcf/1kG_chr{}_vep_pickorder_pantro6filtered_regfiltered_modified_codon.vcf.gz".format(chr))
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
        print("prework done")

        pre_pos = 0
        pre_conseq = None
        pre_flank = None
        pre_pair = None
        pre_pos = 0
        while True:
            try:
                record = next(f1)
            except OSError:
                continue
            except StopIteration:
                break
            # if record.pos in skip_list:
            #     position_skipped+=1
            #     continue
            try:
                vep = record.info["CSQ"]
                # vep = record.info["vep"]
            except KeyError:
                continue
            cur_pos = record.pos
            if cur_pos-pre_pos<=1:
                pre_pos = cur_pos
                position_skipped +=1
                try:
                    if pre_conseq == con:
                        del super_dict[chr][pre_flank][pre_conseq][pre_pos]
                        continue
                    elif pre_conseq == exp:
                        del super_dict[chr][pre_conseq][pre_conseq][pre_pair][pre_pos]
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
            if conseq == con:
                con_processed +=1
                if conseq in super_dict[chr][flanking]:
                    super_dict[chr][flanking][conseq][cur_pos] = ac
                    pre_flank = flanking
                    pre_conseq = conseq
                    pre_pos = cur_pos
                    continue
                elif conseq not in super_dict[chr][flanking]:
                    super_dict[chr][flanking][conseq] = {}
                    super_dict[chr][flanking][conseq][cur_pos] = ac
                    pre_flank = flanking
                    pre_conseq = conseq
                    pre_pos = cur_pos
                    # print(random.choice(zero_one))
                    continue
            elif conseq == exp:
                cur_pair = annotation.split(sep="|")[16]
                exp_processed+=1
                if conseq in super_dict[chr][flanking]:
                    if cur_pair in super_dict[chr][flanking][conseq]:
                        super_dict[chr][flanking][conseq][cur_pair][cur_pos] = ac
                        pre_flank = flanking
                        pre_conseq = conseq
                        pre_pos = cur_pos
                        pre_pair = cur_pair
                        continue
                    elif cur_pair not in super_dict[chr][flanking][conseq]:
                        super_dict[chr][flanking][conseq][cur_pair] = {}
                        super_dict[chr][flanking][conseq][cur_pair][cur_pos] = ac
                        pre_flank = flanking
                        pre_conseq = conseq
                        pre_pos = cur_pos
                        pre_pair = cur_pair
                        continue
                elif conseq not in super_dict[chr][flanking]:
                    super_dict[chr][flanking][conseq] = {}
                    super_dict[chr][flanking][conseq][cur_pair] = {}
                    super_dict[chr][flanking][conseq][cur_pair][cur_pos] = ac
                    pre_flank = flanking
                    pre_conseq = conseq
                    pre_pos = cur_pos
                    pre_pair = cur_pair
                    # print(random.choice(zero_one))
                    continue
            else:
                conseq_skipped +=1
                continue
        print(f" done, processed con:{con_processed}, exp:{exp_processed}, skipped position:{position_skipped}, format:{format_skipped}, sequence:{sequence_skipped}, key_error:{key_error_skipped}, conseq:{conseq_skipped}")
    pickle.dump(super_dict, open(output_fn, 'wb'))
