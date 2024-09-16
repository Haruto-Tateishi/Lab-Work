# this scripts provides relative cumulative SFSs

import numpy as np
# from matplotlib import colormaps
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import random
import pickle

colors = ['black', 'blue', 'red', 'maroon', 'orange', 'greenyellow', 'green', 'cyan', 'violet', 'tan', 'gold', 'darkgrey', 'deepskyblue', 'darkorchid', 'salmon', 'tab:olive', 'tab:brown', 'mediumseagreen', 'steelblue', 'indigo', 'deeppink', 'teal', 'cornflowerblue', 'khaki', 'peru', 'palevioletred', 'lightgrey', 'rosybrown', 'darkkhaki', 'yellow', 'palegreen', 'aquamarine', 'darkgreen', 'lightcyan', 'lightskyblue', 'lightslategrey', 'mediumvioletred', 'pink', 'darkred', 'chocolate', 'firebrick', 'purple', 'cadetblue', 'mediumaquamarine', 'lime', 'lemonchiffon', 'wheat', 'whitesmoke', 'tomato', 'thistle', 'plum', 'royalblue', 'lavender', 'mediumturquoise']
def curve_color(colors):
    color = colors[0]
    colors.remove(color)
    # print(colors)
    return color


def rel_cum(conseq, data_line, curve_color):
    data = data_line.strip().split()
    float_list = [float(x) for x in data]
    np_array = np.array(float_list)
    # print(np_array.dtype)
    cumulative_sum = np.cumsum(np_array)
    # print(cumulative_sum)
    relative_cumulative = cumulative_sum / cumulative_sum[-1]
    x_axis = np.arange(1, len(relative_cumulative) + 1)
    plt.plot(x_axis, relative_cumulative, marker='', linestyle='-', color=curve_color, label=f'{conseq}')

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


# fn1 = "Document/SFS/1kG-vep-conseq-single-donwsample-200.txt"
# f1 = open(fn1, 'r')
fn1 = "Document/Pickle/1kG_sift4g_single_conseq.p"
nc = 5008

# exception_list = ["mature_miRNA_variant", "non_coding_transcript_variant", "coding_transcript_variant", "transcript_ablation", "frameshift_variant", "inframe_insertion", 'stop_lost', "start_lost", 'splice_donor_variant', "stop_retained_variant", "inframe_deletion", "transcript_amplification", "splice_acceptor_variant", "stop_gained"]
exception_list = []

# num_curves = 25
# colormap = get_cmap('gist_ncar', num_curves)

plt.figure(figsize=(10, 6))
plt.xlabel('Data')
plt.ylabel('Relative Cumulative Sum')
plt.title(f'Relative Cumulative Chart with nc value = {str(nc)}')
plt.grid(True)

# lines = f1.readlines()
# input_len = len(lines)
# i = 0
# g = 0
# while True:
#     if i == input_len:
#         break 
#     else:
#         conseq = lines[i].strip()
#         # print(conseq)
#         # if "&" in conseq:
#         #     i +=2
#         #     continue
#         i +=1
#         data_line = lines[i]
#         # print(data_line)
#         i +=2
#         color = curve_color(colors)
#         # color = colormap(g)
#         g += 1
#         if conseq in exception_list:
#             continue
#         rel_cum(conseq, data_line, color)
#         pass


pkl_rel_cum(fn1, nc, colors, exception_list)

plt.ylim(0.85, 1.0001)
plt.legend(loc='lower right', fontsize='small', ncol=3, frameon=True)
plt.show()