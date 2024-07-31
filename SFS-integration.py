from PRF_Ratios import PRF_Ratios_functions as prf
import numpy as np
import scipy.integrate as integrate
import scipy.special as special
import pickle 

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
    data_list = data_line.strip().split()
    float_list = [float(x) for x in data_list]
    np_array = np.array(float_list)
    sum = np.cumsum(np_array)
    try:
        relative = sum / sum[-1]
        result = integrate.cumulative_trapezoid(relative, x=None, dx=1.0, axis=-1, initial=None)
        auc = result[-1]

        line = [conseq, str(auc), "\n"]
        line = "\t".join(line)
        output_file.write(line)
    except RuntimeWarning:
        line = [conseq, "No Curve", "\n"]
        line = "\t".join(line)
        output_file.write(line)

def dict_integration(pickle_file, nc, output_file):
    Dic = pickle.load(open(pickle_file, "rb"))
    for conseq in Dic:
        sfs = np.zeros(nc, dtype=np.float32)
        for acan in Dic[conseq]:
            # an = acan.split(sep="_")[3]
            ac = int(acan.split(sep="_")[1])
            sfs[ac-1] += Dic[conseq][acan]
        np_array = np.array(sfs)
        sum = np.cumsum(np_array)
        try:
            relative = sum / sum[-1]
            result = integrate.cumulative_trapezoid(relative, x=None, dx=1.0, axis=-1, initial=None)
            auc = result[-1]

            line = [conseq, str(auc), "\n"]
            line = "\t".join(line)
            output_file.write(line)
        except RuntimeWarning:
            line = [conseq, "No Curve", "\n"]
            line = "\t".join(line)
            output_file.write(line)



fn1 = "gnomad-SFS-downsample-integrate-single-1000000.tsv"
f1 = open(fn1, "w")
f1.write("conseq\tarea-under-curve\n")

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
        if "&" in conseq:
            i +=2
            continue
        i +=1
        data_line = lines[i]
        i +=1
        data_integration(conseq, data_line, f1)
        pass

f1.close()
f2.close()


# fn1 = "1kG-SFS-syn-single-integrate.tsv"
# f1 = open(fn1, "w")
# f1.write("conseq\tarea-under-curve\n")

# fn2 = "1kG_vep_syn_single_codon_counts.p"

# dict_integration(fn2, 5008, f1)

# f1.close()