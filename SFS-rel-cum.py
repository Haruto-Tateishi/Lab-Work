# this scripts provides relative cumulative SFSs

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def curve_color():
    # pick up color that does not overlap
    return None

def rel_cum(conseq, data_line, curve_color):
    data = data_line.strip().split()
    float_list = [float(x) for x in data]
    np_array = np.array(float_list)
    # print(np_array.dtype)
    cumulative_sum = np.cumsum(np_array)
    # print(cumulative_sum)
    relative_cumulative = cumulative_sum / cumulative_sum[-1]
    x_axis = np.arange(1, len(relative_cumulative) + 1)
    plt.plot(x_axis, relative_cumulative, marker='', linestyle='-', color=f'{curve_color}', label=f'{conseq}')

fn1 = "gnomad-vep-downsample-1000000.txt"
f1 = open(fn1, "r")

lines = f1.readlines()
lines_len = len(lines)
nc = lines[0].strip()

plt.figure(figsize=(10, 6))
plt.xlabel('Data')
plt.ylabel('Relative Cumulative Sum')
plt.title(f'Relative Cumulative Chart with nc value = {str(nc)}')
plt.grid(True)

i = 2
while True:
    if i == lines_len:
        break
    else:
        conseq = lines[i]
        i +=1
        data_line = lines[i]
        i +=1
        color = curve_color()
        rel_cum(conseq, data_line, color)
        pass


plt.legend(loc='lower right', fontsize='small', ncol=3, frameon=True)
plt.show()