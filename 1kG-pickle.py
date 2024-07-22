# this script reads the binary gzipped file, creats a super dictionary where the first keys are consequences, the second keys are AC and AN, and values are the frequency, and saves the dictionary to pickle file.

import gzip
import pysam
import pickle
import os.path as op


super_dict = {}

for i in range(1, 2):
    print("start ",i,end="")
    fn1 = op.join("/home/haruto/ensembl-vep/output.vcf","1kG-chr{}-vep-every.vcf.gz".format(i))
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
        for annotation in vep:
            conseq = annotation.split(sep="|")[1]
            if "&" in annotation: # skip all records that has consequences with "&"
                break
            ac = record.info["AC"][0]
            an = record.info["AN"]
            data = "AC_" + str(ac) + "_AN_" + str(an)
            if conseq not in super_dict:
                super_dict[conseq] = {}
                super_dict[conseq][data] = 1
            else:
                if data not in super_dict[conseq]:
                    super_dict[conseq][data] = 1
                else:
                    super_dict[conseq][data] += 1
    print(" done")

fn2 = op.join("/home/haruto/Lab-Work","1kG_vep_consequence_ANAC_counts.p")
# fn2 = op.join("/mnt/d/genemod/better_dNdS_models/popgen/Human/gnomad/pickles","temp.p")
pickle.dump(super_dict, open(fn2, 'wb'))