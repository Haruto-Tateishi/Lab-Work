# this script is used to downsample the pickle file and generate the SFS

import os.path as op
import numpy as np
from scipy import special
import sys
from collections import defaultdict
import os.path as op 
import pickle
import warnings

# Set up a custom warning handler to catch RuntimeWarnings
def warning_handler(message, category, filename, lineno, file=None, line=None):
    print(f"RuntimeWarning encountered: {message}")
    print(f"In file: {filename}, line {lineno}")
    sys.exit(1)  # Exit the program

warnings.filterwarnings('error', category=RuntimeWarning)
warnings.showwarning = warning_handler

# log_binomial is used to calculate the log of the binomial coefficient
def log_binomial(n, k):
    return special.gammaln(n + 1) - special.gammaln(k + 1) - special.gammaln(n - k + 1)

# hypergeometric is used to calculate the hypergeometric probability
def hypergeometric(N, i, n, k):
    """
        N bins,  for bin i,  probability of sampling k items with sample size n
    """
    try:
        log_p = log_binomial(i, k) + log_binomial(N - i, n - k) - log_binomial(N, n)
    except RuntimeWarning:
        print(f"Exiting due to RuntimeWarning for inputs: N={N}, i={i}, n={n}, k={k}")
        sys.exit(1)        

    if log_p == -np.inf:
        return 0.0
    if log_p > -700:  # exp(-700) is approximately the smallest positive float
        return np.exp(log_p)
    else:
        return 0.0  # Return 0 for extremely small probabilities

# getcodonpairs is used to get the codon pairs and their corresponding amino acids
def getcodonpairs():
    codonpairstrings = ["GCC    GCT","GCA   GCG","GCA   GCC","GCC   GCG","GCA   GCT","GCG   GCT","TGC   TGT","GAC   GAT","GAA   GAG","TTC   TTT","GGC   GGT","GGA   GGC","GGA   GGG","GGA   GGT","GGC   GGG","GGG   GGT","CAC   CAT","ATC   ATT","ATA   ATC","ATA   ATT","AAA   AAG","CTG   TTG","CTA   CTG","CTC   CTT","TTA   TTG","CTG   CTT","CTC   CTG","CTA   TTA","CTA   CTT","CTA   CTC","AAC   AAT","CCA   CCG","CCC   CCT","CCA   CCC","CCC   CCG","CCG   CCT","CCA   CCT","CAA   CAG","CGC   CGT","CGA   CGG","AGA   AGG","CGA   CGC","CGC   CGG","AGG   CGG","AGA   CGA","CGG   CGT","CGA   CGT","TCC   TCT","AGC   AGT","TCA   TCG","TCA   TCC","TCC   TCG","TCG   TCT","TCA   TCT","ACC   ACT","ACA   ACG","ACA   ACC","ACC   ACG","ACA   ACT","ACG   ACT","GTA   GTG","GTC   GTT","GTG   GTT","GTC   GTG","GTA   GTT","GTA   GTC","TAC   TAT"]
    aacodes = ["A","A","A","A","A","A","C","D","E","F","G","G","G","G","G","G","H","I","I","I","K","L","L","L","L","L","L","L","L","L","N","P","P","P","P","P","P","Q","R","R","R","R","R","R","R","R","R","S","S","S","S","S","S","S","T","T","T","T","T","T","V","V","V","V","V","V","Y"]
    return codonpairstrings,aacodes

# build_sampled_SFS is used to build the sampled site frequency spectrum (SFS)
def build_sampled_SFS(consequence,cD,nc):
    """
    dd is a dictionary that has 
    """
    sampledsfs = np.zeros(nc, dtype=np.float32)
    numvals = 0
    numokvals = 0
    skipped_records = 0
    zero_one = [0, 1]
    print(consequence)
    for acanstr in cD:
        temp = acanstr.split(sep="_")
        (ac, an) = (int(temp[1]), int(temp[3]))
        numvals += 1
        if an >= nc:
            # print(random.choice(zero_one))
            numokvals += 1
            popp = ac / an
            for k in range(min(nc, ac + 1)):
                p = hypergeometric(an, ac, nc, k)
                if np.isnan(p):
                    print(f"Warning: NaN result for N={an}, i={ac}, n={nc}, k={k}")
                sampledsfs[k] += cD[acanstr] * p
                if (k / nc > popp) and p < 1e-100:
                    break
        else:
            skipped_records += 1
    print(f"Conseq: {consequence}, Total records: {numvals}, Skipped records: {skipped_records}, Processed records: {numokvals}")
    return sampledsfs, numvals, numokvals

# processnc_dict_pickle is used to process the dictionary and generate the SFS
def processnc_dict_pickle(consequence,cD,nc,sfsoutfilename):
    codonpairlist,aalist = getcodonpairs()
    sumnumvals = 0
    sumnumokvals = 0
    sampledsfs,numvals,numokvals = build_sampled_SFS(consequence,cD,nc)
    sfsf = open(sfsoutfilename, 'a')
    sfsf.write("{}\n".format(consequence))
    sfsf.write(' '.join("{:.1f}".format(x) for x in sampledsfs) + "\n\n")
    sfsf.close()
    return int(sumnumvals/(2*len(codonpairlist))),int(sumnumokvals/(2*len(codonpairlist)))

# processnc_superdict_pickle_all is used to process the super dictionary and generate the SFS
def processnc_superdict_pickle_all(consequence,cD,nc,sfsoutfilename):
    codonpairlist,aalist = getcodonpairs()
    sumnumvals = 0
    sumnumokvals = 0
    sampledsfs,numvals,numokvals = build_sampled_SFS(consequence,cD,nc)
    sfsf = open(sfsoutfilename, 'a')
    sfsf.write("{}\n".format(consequence))
    sfsf.write(' '.join("{:.1f}".format(x) for x in sampledsfs) + "\n\n")
    sfsf.close()
    return int(sumnumvals/(2*len(codonpairlist))),int(sumnumokvals/(2*len(codonpairlist)))

# downsample_dict_pickle is used to downsample the dictionary and generate the SFS
def downsample_dict_pickle():
    nc = 5008
    picklefilename = "Document/Pickle/FlankingBases/syn_codon/1kG_vep_pickorder_flanking_syn_codon.p"
    sfsoutfilename = "uk10k_vep_pickorder_flanking_CpGfiltered_intergenic_synonymous_down_5008.txt"
    # sfsoutfile = open(sfsoutfilename,'w')
    # sfsoutfile.write("{}{}{}{}".format(nc,"\n",picklefilename,"\n"))
    # sfsoutfile.close()
    Dictionary = pickle.load(open(picklefilename, 'rb'))

    for first_key in Dictionary:
        processnc_dict_pickle(first_key,Dictionary[first_key],nc,sfsoutfilename)

# downsample_superdict_pickle is used to downsample the super dictionary and generate the SFS
def downsample_superdict_pickle():
    nc = 30
    picklefilename = "Document/Pickle/FlankingBases/syn_codon/1kG_vep_pickorder_flanking_syn_codon_paired.p"
    sfsoutfilename = f"1kG_vep_pickorder_flanking_syn_codon_all_down_{nc}_2.txt"
    sfsoutfile = open(sfsoutfilename,'w')
    # sfsoutfile.write("{}{}{}{}".format(nc,"\n",picklefilename,"\n"))
    sfsoutfile.close()
    Dictionary = pickle.load(open(picklefilename, 'rb'))
    ind = 0
    for first_key in Dictionary:
        syn_codon = list(first_key)
        syn_codon[3] = "-"
        syn_codon = "".join(syn_codon)
        # ind +=1
        # if ind == 5:
        #     break
        for second_key in Dictionary[first_key]:
            sfsoutfile = open(sfsoutfilename, 'a')
            sfsoutfile.write(f"{first_key}\n")
            sfsoutfile.close()
            processnc_superdict_pickle_all(second_key,Dictionary[first_key][second_key],nc,sfsoutfilename)

# Main function to run the downsampling
downsample_superdict_pickle()