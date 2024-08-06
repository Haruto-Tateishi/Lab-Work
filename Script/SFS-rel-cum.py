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

def pkl_rel_cum(pickle_file, nc, colormap):
    Dic = pickle.load(open(pickle_file, "rb"))
    i = 0
    for conseq in Dic:
        sfs = np.zeros(nc, dtype=np.float32)
        c_color = colormap(i)
        # print(c_color)
        i +=1
        for acan in Dic[conseq]:
            # an = acan.split(sep="_")[3]
            ac = int(acan.split(sep="_")[1])
            sfs[ac-1] += Dic[conseq][acan]
        np_array = np.array(sfs)
        sum = np.cumsum(np_array)
        try:
            relative = sum / sum[-1]
            x_axis = np.arange(1, len(relative) + 1)
            plt.plot(x_axis, relative, marker='', linestyle='-', color=c_color, label=f'{conseq}')
        except RuntimeWarning:
            pass


# fn1 = "1kG-vep-syn-single-codon-downsample-500.txt"
# f1 = open(fn1, 'r')
fn1 = "1kG_vep_syn_all_codon_counts.p"
nc = 5008

num_curves = 134
colormap = get_cmap('gist_ncar', num_curves)

plt.figure(figsize=(10, 6))
plt.xlabel('Data')
plt.ylabel('Relative Cumulative Sum')
plt.title(f'Relative Cumulative Chart with nc value = {str(nc)}')
plt.grid(True)

# lines = f1.readlines()
# input_len = len(lines)
# i = 2
# g = 0
# while True:
#     if i == input_len:
#         break 
#     else:
#         conseq = lines[i].strip()
#         # if "&" in conseq:
#         #     i +=2
#         #     continue
#         i +=1
#         data_line = lines[i]
#         i +=1
#         # color = curve_color(colors)
#         color = colormap(g)
#         g += 1
#         rel_cum(conseq, data_line, color)
#         pass

pkl_rel_cum(fn1, 5008, colormap)

plt.ylim(0.9, 1.0001)
plt.legend(loc='lower right', fontsize='small', ncol=8, frameon=True)
plt.show()