# this script reads the binary gzipped file, creats a super dictionary where the first keys are consequences, the second keys are AC and AN, and values are the frequency, and saves the dictionary to pickle file.

import gzip
import pysam
import pickle

fn1 = "toy-1kG-chr22-vep-pick.vcf.bin.gz"
gunzipped = gzip.open(fn1, 'rt') 
f1 = pysam.VariantFile(gunzipped)


conseq_list = []
super_dict = {}

for record in f1:
    # print(record)
    try:
        vep = record.info["CSQ"][0]
    except KeyError:
        continue
    conseq = vep.split(sep="|")[1]
    ac = record.info["AC"][0]
    an = record.info["AN"]
    if conseq not in conseq_list:
        conseq_list.append(conseq)
        data = "AC_" + str(ac) + "_AN_" + str(an)
        super_dict[conseq] = {}
        super_dict[conseq][data] = 1
        continue
    if conseq in conseq_list:
        data = "AC_" + str(ac) + "_AN_" + str(an)
        dict_keys = super_dict[conseq].keys()
        if data not in dict_keys:
            super_dict[conseq][data] = 1
            continue
        if data in dict_keys:
            super_dict[conseq][data] += 1
            continue

# print(super_dict)

fn2 = "toy-1kG-chr22-dict.pkl"
pickle.dump(super_dict, open(fn2, 'wb'))