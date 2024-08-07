# this script simulates SFSs with different parameters for simsfs() function.

from PRF_Ratios import PRF_Ratios_functions as prf
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

a, b = prf.simsfs(20, 100, 1000, None, 900, True)
print(a, b)

np_array_a = np.array(a)[1:]
np_array_b = np.array(b)[1:]

sum_a = np.cumsum(np_array_a)
sum_b = np.cumsum(np_array_b)

relative_a = sum_a / sum_a[-1]
relative_b = sum_b / sum_b[-1]
# print(relative_x, relative_y)

x_axis_a = np.arange(1, len(relative_a) + 1)
x_axis_b = np.arange(1, len(relative_b) + 1)


plt.figure(figsize=(10, 6))
plt.plot(x_axis_a, relative_a, marker='', linestyle='-', color='b', label='Data')
plt.xlabel('Data')
plt.ylabel('Relative Cumulative Sum')
plt.title('Relative Cumulative Chart')
plt.grid(True)
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(x_axis_b, relative_b, marker='', linestyle='-', color='b', label='Data')
plt.xlabel('Data')
plt.ylabel('Relative Cumulative Sum')
plt.title('Relative Cumulative Folded Chart')
plt.grid(True)
plt.legend()
plt.show()