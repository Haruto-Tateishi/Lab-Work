# This script reads the binary gzipped file, creats a super dictionary.
# Due to lack of memory in my computer, the super dictionary is saved in chunks.

import gzip
import pysam
import pickle
import os
import random
import paramiko
import io
import logging
import gc

# Set up logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# list all files in the directory
def list_documents(directory, extension):
    try:
        file_names = os.listdir(directory)
        file_names = sorted(file_names)
        document_extensions = {extension}
        documents = [f for f in file_names if os.path.isfile(os.path.join(directory, f)) and os.path.splitext(f)[1] in document_extensions]
        return documents
    except Exception as e:
        logging.error(f"Error listing documents: {e}", exc_info=True)

# list all files in the directory with sftp
def list_directory(sftp, directory):
    """List the contents of a directory on the remote server."""
    try:
        return sftp.listdir(directory)
    except Exception as e:
        print(f"Error listing directory: {e}")
        return []

# read the file in chunks
def read_remote_file_in_chunks(sftp_file, chunk_size=1024):
    """Generator to read a remote file in chunks."""
    while True:
        data = sftp_file.read(chunk_size)
        if not data:
            break
        yield data


# read the gzipped file line by line
def read_gzipped_file_line_by_line(file_path):
    with gzip.open(file_path, 'rt') as gz_file:
        for line in gz_file:
            yield line.strip()

# generate a super dictionary from roulette vcf file which contains mutation rate info.
def mutation_rate_pickle(input_fn, super_dict, chr, bank):
    zero_one = [0,1]

    print("program start")
    total_lines_processed = 0
    chunk_counter = 0
    chunk_size = 100000000
    for line in read_gzipped_file_line_by_line(input_fn):
        if len(line) <1:
            break
        # print(random.choice(zero_one))
        # print(record.info["MR"])
        if line[0] =="#":
            continue
        line = line.split(sep="\t")
        mr = line[7].split(sep=";")[1]
        rate = mr.split(sep="=")[1]
        pos = line[1]
        mut = f"{line[3]}-{line[4]}"
        # print(pos, mut, rate)
        if pos in super_dict.keys(): 
            super_dict[pos][mut] = rate
        else:
            super_dict.update({pos:{mut:rate}})
            # print(rate)
        
        total_lines_processed += 1
        # print(total_lines_processed)

        if total_lines_processed % chunk_size == 0:
            chunk_counter += 1
            # Save intermediate super_dict to pickle file
            output_fn = f"{bank}_chr{chr}_mutation_rate_sorted_{chunk_counter}.p"
            with open(output_fn, 'wb') as f:
                pickle.dump(super_dict, f)
                f.close()
            # Clear memory
            del super_dict
            gc.collect()
            super_dict = {}
            total_lines_processed = 0
        

    # Save any remaining data
    if super_dict:
        chunk_counter += 1
        output_fn = f"{bank}_chr{chr}_mutation_rate_sorted_{chunk_counter}.p"
        with open(output_fn, 'wb') as f:
            pickle.dump(super_dict, f)
            f.close()
        del super_dict
        gc.collect()
        super_dict = {}
        total_lines_processed = 0

# assign GC bias based on the reference and alternate alleles
def GC_bias_test(ref, alt):
    gc = ["G", "C"]
    at = ["A", "T"]
    if ref in gc and alt in gc or ref in at and alt in at:
        return "n"
    elif ref in at and alt in gc:
        return "A"
    elif ref in gc and alt in at:
        return "B"
    
# fetch the bin number pickle file 
def bin_num_pickle_fetch(std_chr):
    bin_num_directory = "Document/Pickle/MutationRate/SuperDictionary/BinNumber"
    documents = list_documents(bin_num_directory, extension=".p")
    for bin_num_pickle in documents:
        if bin_num_pickle.split(sep="_")[1] != std_chr:
            continue
        else:
            bin_num_pickle = os.path.join(bin_num_directory, bin_num_pickle)
            with open(bin_num_pickle, "rb") as f:
                bin_num_dict = pickle.load(f)
            return bin_num_dict

# synonymous variant
#     GC bias
#         codon pair
#             mutation rate bin number
#                 position
#                     an_ac
def GC_bias_synonymous(vcf, std_conseq, std_chr, bin_num_pickle):
    zero_one = [0, 1]
    super_dict = {}
    super_dict.update({"A":{}, "B":{}, "n":{}})
    for record in vcf:
        vep = record.info["CSQ"][0]
        conseq = vep.split(sep="|")[1]
        if conseq != std_conseq:
            # print(conseq)
            continue
        pos = record.pos
        an_ac = str(record.info["AN"]) + "_" + str(record.info["AC"][0])
        ref = record.ref
        alt = record.alts[0]
        mut = f"{ref}-{alt}"
        GC_bias = GC_bias_test(ref, alt)
        codon_pair = vep.split(sep="|")[16]        
        # print(f"pos:{pos}:{type(pos)}, mut:{mut}:{type(mut)}, ref:{ref}:{type(ref)}, alt:{alt}:{type(alt)}, an_ac:{an_ac}:{type(an_ac)}, ")
        try:
            bin_num = bin_num_pickle[str(pos)][mut]
            if bin_num == 14:
                continue
        except KeyError:
            continue
        # generate superdictionary from here
        if codon_pair in super_dict[GC_bias].keys():
            if bin_num in super_dict[GC_bias][codon_pair]:
                super_dict[GC_bias][codon_pair][bin_num][pos] = an_ac
            else:
                super_dict[GC_bias][codon_pair].update({bin_num:{pos:an_ac}})
        else:
            super_dict[GC_bias].update({codon_pair:{bin_num:{pos:an_ac}}})
    return super_dict

# intergenic variant or missense
#     GC bias
#         mutation rate bin number
#             position
#                 an_ac
def GC_bias_intergenic_missense(vcf, std_conseq, std_chr, bin_num_pickle):
    zero_one = [0, 1]
    super_dict = {}
    super_dict.update({"A":{}, "B":{}, "n":{}})
    for record in vcf:
        vep = record.info["CSQ"][0]
        conseq = vep.split(sep="|")[1]
        if conseq != std_conseq:
            # print(conseq)
            continue
        pos = record.pos
        an_ac = str(record.info["AN"]) + "_" + str(record.info["AC"][0])
        ref = record.ref
        alt = record.alts[0]
        mut = f"{ref}-{alt}"
        GC_bias = GC_bias_test(ref, alt)
        # print(f"pos:{pos}:{type(pos)}, mut:{mut}:{type(mut)}, ref:{ref}:{type(ref)}, alt:{alt}:{type(alt)}, an_ac:{an_ac}:{type(an_ac)}, ")
        try:
            bin_num = bin_num_pickle[str(pos)][mut]
            if bin_num == 14:
                continue
        except KeyError:
            continue
        # generate superdictionary from here
        if bin_num in super_dict[GC_bias].keys():
            super_dict[GC_bias][bin_num][pos] = an_ac
        else:
            super_dict[GC_bias].update({bin_num:{pos:an_ac}})
    return super_dict

# generate a super dictionary about GC bias from the vcf file
def super_dict_GC_bias(input_fn, std_conseq_list, std_chr, bin_num_pickle, bank):
    for std_conseq in std_conseq_list:
        output_fn = f"{bank}_{std_chr}_CG_{std_conseq[0:3]}.p"
        with open(input_fn, 'rb') as f:
            vcf = pysam.VariantFile(f)
        print("VCF is fetched")
        if std_conseq == "synonymous_variant":
            print(f"dealing with {std_conseq}")
            super_dict = GC_bias_synonymous(vcf, std_conseq, std_chr, bin_num_pickle)
        elif std_conseq == "intergenic_variant" or std_conseq == "missense_variant":
            print(f"dealing with {std_conseq}")
            super_dict = GC_bias_intergenic_missense(vcf, std_conseq, std_chr, bin_num_pickle)
        pickle_save(output_fn, super_dict)

# save the super dictionary to pickle file 
def pickle_save(output_fn, obj):
    with open(output_fn, "wb") as f:
        pickle.dump(obj, f)

# main function for GC bias super dictionary
def GC_bias_script():
    std_conseq_list = ["synonymous_variant", "intergenic_variant", "missense_variant"]
    # Select gnomad, uk10k or 1kG
    bank = "uk10k"
    directory = 'vcf'
    documents = list_documents(directory, extension=".bgz")
    for chr in range(21,22):
        std_chr = f"chr{chr}"
        print(f"{std_chr} started")
        bin_num_pickle = bin_num_pickle_fetch(std_chr)
        print("Dictionary is fetched")
        for input_fn in documents:
            input_fn = os.path.join(directory, input_fn)
            if bank == "gnomad":
                cmp_chr = input_fn.split(sep=".")[5]
            elif bank == "uk10k":
                cmp_chr = input_fn.split(sep="_")[1]
            elif bank == "1kG":
                cmp_chr = input_fn.split(sep=".")[1]
            if cmp_chr != std_chr:
                continue
            super_dict_GC_bias(input_fn, std_conseq_list, std_chr, bin_num_pickle, bank)



GC_bias_script()