# this scripts provides relative cumulative SFSs

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


fn1 = "debugdownsampling-1000000.txt"
f1 = open(fn1, "r")

lines = f1.readlines()
data = lines[3].strip().split()

float_list = [float(x) for x in data]
np_array = np.array(float_list)
print(np_array.dtype)

cumulative_sum = np.cumsum(np_array)
print(cumulative_sum)

relative_cumulative = cumulative_sum / cumulative_sum[-1]

x_axis = np.arange(1, len(relative_cumulative) + 1)

plt.figure(figsize=(10, 6))
plt.plot(x_axis, relative_cumulative, marker='', linestyle='-', color='b', label='Data')
plt.xlabel('Data')
plt.ylabel('Relative Cumulative Sum')
plt.title('Relative Cumulative Chart')
plt.grid(True)
plt.legend()
plt.show()