# this script simulates SFSs with different parameters for simsfs() function.

from PRF_Ratios import PRF_Ratios_functions as prf
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

a, b = prf.simsfs(0.1, 1, 100000, None, 90000, False)
print(a, b)

np_array_a = np.array(a)
np_array_b = np.array(b)

sum_a = np.cumsum(np_array_a)
sum_b = np.cumsum(np_array_b)

relative_a = sum_a / sum_a[-1]
relative_b = sum_b / sum_b[-1]
# print(relative_x, relative_y)

x_axis_a = np.arange(1, len(relative_a) + 1)
x_axis_b = np.arange(1, len(relative_b) + 1)

def logarithmic_fit(x, a, b):
    return a * np.log(x) + b

params_a, params_covariance_a = curve_fit(logarithmic_fit, x_axis_a, relative_a)
curve_a = logarithmic_fit(x_axis_a, *params_a)

params_b, params_covariance_b = curve_fit(logarithmic_fit, x_axis_b, relative_b)
curve_b = logarithmic_fit(x_axis_b, *params_b)


plt.figure(figsize=(10, 6))
plt.plot(x_axis_a, relative_a, marker='', linestyle='-', color='b', label='Data')
plt.plot(x_axis_a, curve_a, linestyle='-', color='r', label='Curve Fit')
plt.xlabel('Data')
plt.ylabel('Relative Cumulative Sum')
plt.title('Relative Cumulative Chart')
plt.grid(True)
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(x_axis_b, relative_b, marker='', linestyle='-', color='b', label='Data')
plt.plot(x_axis_b, curve_b, linestyle='-', color='r', label='Curve Fit')
plt.xlabel('Data')
plt.ylabel('Relative Cumulative Sum')
plt.title('Relative Cumulative Folded Chart')
plt.grid(True)
plt.legend()
plt.show()