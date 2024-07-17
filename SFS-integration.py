from PRF_Ratios import PRF_Ratios_functions as prf
import numpy as np
import scipy.integrate as integrate
import scipy.special as special

def sim_integration(theta, g, nc, misspec, maxi, returnexpected, output_file):
    a_list = []
    b_list = []
    for i in range(10):
        a, b = prf.simsfs(theta, g, nc, misspec, maxi, returnexpected)
        np_array_a = np.array(a)
        np_array_b = np.array(b)

        sum_a = np.cumsum(np_array_a)
        sum_b = np.cumsum(np_array_b)

        try:
            relative_a = sum_a / sum_a[-1]
            relative_b = sum_b / sum_b[-1]
        except RuntimeWarning:
            continue

        result_a = integrate.cumulative_trapezoid(relative_a, x=None, dx=1.0, axis=-1, initial=None)
        result_b = integrate.cumulative_trapezoid(relative_b, x=None, dx=1.0, axis=-1, initial=None)

        integrate_a = result_a[-1]
        integrate_b = result_b[-1]

        a_list.append(integrate_a)
        b_list.append(integrate_b)

    average_a = np.average(a_list)
    average_b = np.average(b_list)
    print(average_a, average_b)

    line = [str(theta), str(g), str(nc), str(misspec), str(maxi), str(returnexpected), str(average_a), str(average_b), "\n"]
    line = "\t".join(line)
    output_file.write(line)

def data_integration(conseq, data_line, output_file):
    data_list = data_line.strip().split(sep="\t")
    np_array = np.array(data_list)
    sum = np.cumsum(np_array)
    relative = sum / sum[-1]
    result = integrate.cumulative_trapezoid(relative, x=None, dx=1.0, axis=-1, initial=None)
    integrate = result[-1]

    line = [conseq, str(integrate), "\n"]
    line = "\t".join(line)
    output_file.write(line)

fn1 = "SFS-downsample-integrate.tsv"
f1 = open(fn1, "a")

fn2 = "gnomad-vep-downsample-1000000.txt"
f2 = open(fn2, "r")

lines = f2.readlines()
input_len = len(lines)
i = 2
while True:
    if i == input_len:
        break 
    else:
        conseq = lines[i].strip()
        i +=1
        data_line = lines[i]
        i +=1
        data_integration(conseq, data_line, f1)
        pass

f1.close()
f2.close()