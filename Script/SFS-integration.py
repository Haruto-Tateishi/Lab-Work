from PRF_Ratios import PRF_Ratios_functions as prf
import numpy as np
import scipy.integrate as integrate
import scipy.special as special
import pickle 

def order_check(list_in):
    list_len = len(list_in)
    for i in range(list_len-1):
        # print(list_in[i])
        area_1 = list_in[i][1]
        area_2 = list_in[i+1][1]
        # print(area_1, area_2)
        if area_1 > area_2:
            continue
        if area_2 > area_1:
            return True
    return False

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

def data_integration(input_file, output_file):
    lines = input_file.readlines()
    input_len = len(lines)
    area_list = []
    no_curve_list = []
    i = 0
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
            i +=2
            data_list = data_line.strip().split()
            float_list = [float(x) for x in data_list]
            np_array = np.array(float_list)
            sum = np.cumsum(np_array)
            try:
                relative = sum / sum[-1]
                result = integrate.cumulative_trapezoid(relative, x=None, dx=1.0, axis=-1, initial=None)
                auc = result[-1]
                line = [conseq, auc, "\n"]
                area_list.append(line)
            except RuntimeWarning:
                line = [conseq, "No Curve", "\n"]
                line = "\t".join(line)
                no_curve_list.append(line)

    while True:
        if order_check(area_list):
            for i in range(len(area_list)):
                line_1 = area_list[i]
                area_1 = line_1[1]
                current_pos = i
                for add_pos in range(1, len(area_list)-1):
                    try:
                        line_2 = area_list[i + add_pos]
                        area_2 = line_2[1]
                        if area_1 < area_2:
                            area_list[current_pos] = line_2
                            area_list[i + add_pos] = line_1
                            current_pos = i + add_pos
                            pass
                    except IndexError:
                        break
        else:
           break
   
    for line in area_list:
        line[1] = str(line[1])
        line = "\t".join(line)
        output_file.write(line)


    for line in no_curve_list:
        output_file.write(line)


def dict_integration(pickle_file, nc, output_file, exception_list):
    Dic = pickle.load(open(pickle_file, "rb"))
    area_list = []
    no_curve_list = []
    for conseq in Dic:
        # if "&" in conseq or conseq in exception_list:
        #     continue
        sfs = np.zeros(nc+1, dtype=np.float32)
        for acan in Dic[conseq]:
            # an = acan.split(sep="_")[3]
            ac = int(acan.split(sep="_")[1])
            sfs[ac] += Dic[conseq][acan]
        np_array = np.array(sfs)
        sum = np.cumsum(np_array)
        try:
            relative = sum / sum[-1]
            result = integrate.cumulative_trapezoid(relative, x=None, dx=1.0, axis=-1, initial=None)
            auc = result[-1]
            line = [conseq, auc, "\n"]
            area_list.append(line)
        except RuntimeWarning:
            line = [conseq, "No Curve", "\n"]
            line = "\t".join(line)
            no_curve_list.append(line)
    
    while True:
        if order_check(area_list):
            for i in range(len(area_list)):
                line_1 = area_list[i]
                area_1 = line_1[1]
                current_pos = i
                for add_pos in range(1, len(area_list)-1):
                    try:
                        line_2 = area_list[i + add_pos]
                        area_2 = line_2[1]
                        if area_1 < area_2:
                            area_list[current_pos] = line_2
                            area_list[i + add_pos] = line_1
                            current_pos = i + add_pos
                            pass
                    except IndexError:
                        break
        else:
            break
    
    for line in area_list:
        line[1] = str(line[1])
        line = "\t".join(line)
        output_file.write(line)

    for line in no_curve_list:
        output_file.write(line)



fn1 = "1kG-SFS-vep-syn-mis-all-integrate.tsv"
f1 = open(fn1, "w")
f1.write("conseq\tarea-under-curve\n")

fn2 = "Document/SFS/1kG-vep-syn-mis-all.txt"
f2 = open(fn2, "r")

data_integration(f2, f1)

f1.close()
f2.close()


# fn1 = "1kG-SFS-vep-syn-mis-all-integrate.tsv"
# f1 = open(fn1, "w")
# f1.write("conseq\tarea-under-curve\n")

# fn2 = "/Users/harutotateishi/Library/CloudStorage/OneDrive-TempleUniversity/Documents/Study/Lab/Lab-Work/Document/pickle/1kG_vep_syn_mis_ANAC_counts_all.p"

# exception_list = ["mature_miRNA_variant", "non_coding_transcript_variant", "coding_transcript_variant", "transcript_ablation", "frameshift_variant", "inframe_insertion", 'stop_lost', "start_lost", 'splice_donor_variant', "stop_retained_variant", "inframe_deletion", "transcript_amplification", "splice_acceptor_variant"]
# exception_list = []

# dict_integration(fn2, 5008, f1, exception_list)

# f1.close()