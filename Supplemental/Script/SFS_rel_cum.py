# this scripts provides relative cumulative SFSs

import numpy as np
# from matplotlib import colormaps
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import random
import pickle
import seaborn as sns
import os
import logging

colors = ['black', 'blue', 'red', 'maroon', 'orange', 'greenyellow', 'green', 'cyan', 'violet', 'tan', 'gold', 'darkgrey', 'deepskyblue', 'darkorchid', 'salmon', 'tab:olive', 'tab:brown', 'mediumseagreen', 'steelblue', 'indigo', 'deeppink', 'teal', 'cornflowerblue', 'khaki', 'peru', 'palevioletred', 'lightgrey', 'rosybrown', 'darkkhaki', 'yellow', 'palegreen', 'aquamarine', 'darkgreen', 'lightcyan', 'lightskyblue', 'lightslategrey', 'mediumvioletred', 'pink', 'darkred', 'chocolate', 'firebrick', 'purple', 'cadetblue', 'mediumaquamarine', 'lime', 'lemonchiffon', 'wheat', 'whitesmoke', 'tomato', 'thistle', 'plum', 'royalblue', 'lavender', 'mediumturquoise']
# colors = ['maroon', 'blue', 'black', 'violet', 'green',  'greenyellow', 'orange', 'darkorchid', 'red', 'deepskyblue', 'tan', 'tan', 'cyan', 'tan', 'tan',  'tan', 'tan', 'tan', 'tan', 'tan', 'tan', 'tan', 'tan', 'khaki', 'peru', 'palevioletred', 'lightgrey', 'rosybrown', 'darkkhaki', 'yellow', 'palegreen', 'aquamarine', 'darkgreen', 'lightcyan', 'lightskyblue', 'lightslategrey', 'mediumvioletred', 'pink', 'darkred', 'chocolate', 'firebrick', 'purple', 'cadetblue', 'mediumaquamarine', 'lime', 'lemonchiffon', 'wheat', 'tan', 'tan', 'tan', 'tan', 'tan', 'tan', 'tan']
def curve_color(colors):
    color = colors[0]
    colors.remove(color)
    # print(colors)
    return color

def list_documents(directory, extension):
    try:
        file_names = os.listdir(directory)
        file_names = sorted(file_names)
        document_extensions = {extension}
        documents = [f for f in file_names if os.path.isfile(os.path.join(directory, f)) and os.path.splitext(f)[1] in document_extensions]
        return documents
    except Exception as e:
        logging.error(f"Error listing documents: {e}", exc_info=True)

def rel_cum(input_fn, exception_list, colors, nc):
    fn = open(input_fn, "r")
    lines = fn.readlines()
    len_lines = len(lines)
    i = 0
    g = 0
    while True:
        if i >= len_lines:
            break 
        # if i >= 20:
        #     break
        else:
            # conseq = lines[i].strip() + "_" + lines[i+1].split(sep="_")[0][0:3]
            # syn_codon = lines[i].strip()
            # syn_codon_filename = list(syn_codon)
            # syn_codon_filename[3] = "-"
            # syn_codon_filename = "".join(syn_codon_filename)
            # print(conseq)
            conseq = lines[i]
            if conseq in exception_list:
                color = colors[g]
                i += 3
                g += 1
                continue
            if "&" in conseq:
                i += 3
                continue
            i +=1
            data_line = lines[i]
            # print(data_line)
            i +=1
            if i%4==0:
                color = colors[0]
            else:
                color = colors[1]
            # color = colormap(g)
            # g += 1
            data = data_line.strip().split()
            float_list = [float(x) for x in data]
            np_array = np.array(float_list[1:])
            # print(np_array.dtype)
            cumulative_sum = np.cumsum(np_array)
            # print(cumulative_sum)
            relative_cumulative = cumulative_sum / cumulative_sum[-1]
            x_axis = np.arange(1, len(relative_cumulative) + 1)
            if i%4 ==0:
                plt.plot(x_axis, relative_cumulative, marker='', linestyle='-', color=color, label=f'{conseq}')
                plt.legend(loc='lower right', fontsize=15, ncol=1, frameon=True)
                codon_pair = input_fn.split(sep="_")[3].split(sep=".")[0]
                path = f"Charts/MutationRates/GCBias/1kG_GC_type_{codon_pair}.png"
                # path = f"Charts/FlankingBases/syn_codon/1kG/Cum_1kG_vep_pickorder_flanking_single_{syn_codon_filename}_down_{nc}.png"
                plt.savefig(path)
                plt.close()
                pass
            else:
                plt.figure(figsize=(10, 6))
                plt.xlabel('Allele Counts')
                plt.ylabel('Variant Counts')
                plt.title(f'Relative Cumulative SFS with NC Value = {str(nc)}')
                plt.grid(True)
                plt.ylim(0.9, 1.001)
                plt.plot(x_axis, relative_cumulative, marker='', linestyle='-', color=color, label=f'{conseq}')
                # print("continuing")
                pass

def sfs(input_fn, exception_list, colors):
    fn = open(input_fn, "r")
    lines = fn.readlines()
    len_lines = len(lines)
    i = 0
    g = 0
    while True:
        if i == len_lines:
            break 
        else:
            conseq = lines[i].strip()
            # print(conseq)
            if conseq in exception_list:
                color = colors[g]
                i += 3
                g += 1
                continue
            if "&" in conseq:
                i += 3
                continue
            i +=1
            data_line = lines[i]
            # print(data_line)
            i +=2
            color = colors[g]
            # color = colormap(g)
            g += 1
            data = data_line.strip().split()
            float_list = [float(x) for x in data]
            np_array = np.array(float_list)
            x_axis = np.arange(1, len(np_array) + 1)
            plt.plot(x_axis, np_array, marker='', linestyle='-', color=color, label=f'{conseq}')
            pass

def pkl_rel_cum(pickle_file, nc, colormap, exception_list):
    Dic = pickle.load(open(pickle_file, "rb"))
    i = 0
    for conseq in Dic:
        if "&" in conseq:
            continue
        c_color = colormap[i]
        # print(c_color)
        i +=1
        if conseq in exception_list:
            continue
        sfs = np.zeros(nc+1, dtype=np.float32)
        for acan in Dic[conseq]:
            # an = acan.split(sep="_")[3]
            ac = int(acan.split(sep="_")[1])
            sfs[ac] += Dic[conseq][acan]
        np_array = np.array(sfs)
        sum = np.cumsum(np_array)
        try:
            relative = sum / sum[-1]
            x_axis = np.arange(1, len(relative) + 1)
            plt.plot(x_axis, relative, marker='', linestyle='-', color=c_color, label=f'{conseq}')
        except RuntimeWarning:
            pass

def pkl_sfs(pickle_file, nc, colormap, exception_list):
    Dic = pickle.load(open(pickle_file, "rb"))
    i = 0
    for conseq in Dic:
        if "&" in conseq:
            continue
        c_color = colormap[i]
        # print(c_color)
        i +=1
        if conseq in exception_list:
            continue
        sfs = np.zeros(nc+1, dtype=np.float32)
        for acan in Dic[conseq]:
            # an = acan.split(sep="_")[3]
            ac = int(acan.split(sep="_")[1])
            sfs[ac] += Dic[conseq][acan]
        np_array = np.array(sfs)
        x_axis = np.arange(1, len(np_array) + 1)
        plt.plot(x_axis, np_array, marker='', linestyle='-', color=c_color, label=f'{conseq}')

sfs_path = "Document/SFS/MutationRate/GCBias"
# pkl = "Document/pickle/Flanking_Bases/1kG_vep_pickorder_pantro6filtered_regfiltered_all_integenic_synonymous.p"
nc = 5008
# num_colors = 300
# colors = plt.cm.hsv(np.linspace(0,1,num_colors))
# colors = enumerate(colors)

# exception_list = [ "frameshift_variant", "stop_gained", "regulatory_region_variant", "mature_miRNA_variant", "non_coding_transcript_variant", "coding_transcript_variant", "transcript_ablation",  "inframe_insertion", 'stop_lost', "start_lost", 'splice_donor_variant', "stop_retained_variant", "inframe_deletion", "transcript_amplification", "splice_acceptor_variant", "protein_altering_variant", "coding_sequence_variant"]
exception_list = []


documents = list_documents(sfs_path, ".txt")
# print(documents)
for sfs_fn in documents:
    sfs_fn = os.path.join(sfs_path, sfs_fn)
    rel_cum(sfs_fn, exception_list, colors, nc)

# sfs(sfs_fn, exception_list)

# pkl_rel_cum(pkl, nc, colors, exception_list)

# pkl_sfs(pkl, nc, colors, exception_list)

