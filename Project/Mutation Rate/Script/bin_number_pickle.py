# This script generates the bins of float in which the number of floats in each bin is almost equal

import matplotlib.pyplot as plt
import pickle
import os
import gzip
import heapq
import logging
import pandas as pd

# Set up logging
logging.basicConfig(
    filename='script.log',
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

#  list of all the documents in the directory
def list_documents(directory):
    try:
        file_names = os.listdir(directory)
        file_names = sorted(file_names)
        document_extensions = {".gz"}
        documents = [f for f in file_names if os.path.isfile(os.path.join(directory, f)) and os.path.splitext(f)[1] in document_extensions]
        return documents
    except Exception as e:
        logging.error(f"Error listing documents: {e}", exc_info=True)

# congfigure the graph
def config_graph(num_bins):
    plt.xlabel(f'Mutation Rate {num_bins}bins')
    plt.ylabel('Counts')
    plt.figure(figsize=(10, 6))
    plt.grid(True)
    plt.title(f'Counts of Alleles Based on Mutation Rate')

# create a histogram
def histogram(bins, num_bins, graph_fn):    
    x_axis = [x+1 for x in range(num_bins)]
    plt.plot(x_axis, bins)
    plt.savefig(graph_fn)
    plt.close()

# save the data to a txt file
def data_save_txt(bins, edges, data_fn):
    try:
        bins = [str(x) for x in bins]
        line_bins = " ".join(bins) + "\n"

        # str_edges = [str(y) for y in edges.tolist()]  # Convert NumPy array to list
        # line_edges = " ".join(str_edges) + "\n"

        with open(data_fn, 'w') as f:
            f.write(line_bins)
            # f.write(line_edges)
    except Exception as e:
        logging.error(f"Error in data_save_txt: {e}", exc_info=True)

# save a list of floats to a file
def float_list_save(float_list, data_fn):
    str_float = [str(x) for x in float_list]
    line_float = " ".join(str_float)
    with open(data_fn, "w") as f:
        f.write(line_float)

# Function to insert a float into a min-heap
def insert_into_heap(heap, value):
    heapq.heappush(heap, value)

# Function to convert a min-heap to a sorted list
def heap_to_sorted_list(heap, sorted_list):
    while heap:
        sorted_list.append(heapq.heappop(heap))
    return sorted_list

# print the progress of the dictionary
def print_progress(dict):
    for key in dict:
        print(f"{key}:{dict[key]}")

# Function to create a histogram of mutation rates
def even_histogram_bins(input_fn, mut_dict, float_list):
    try:
        with gzip.open(input_fn, 'rb') as f:
            super_dict = pickle.load(f)
        logging.info("Process Start")
        for pos in super_dict:
            rate_max = 0.0
            for mut in super_dict[pos]:
                if mut not in mut_dict:
                    mut_dict[mut] = 1
                else:
                    mut_dict[mut] +=1 
                rate = float(super_dict[pos][mut])
                float_list.append(rate)
        logging.info("Process End")
        return float_list, mut_dict
    except Exception as e:
        logging.error(f"Error in even_histogram_bins: {e}", exc_info=True)

# Function to create bins based on quantiles
def interval_range(num_bins, float_list):
    try:
        logging.debug(f'Creating pandas Series from float_list with {len(float_list)} elements')
        float_series = pd.Series(float_list)

        logging.debug('Creating quantile-based bins with pd.qcut')
        float_series_binned, bin_edges = pd.qcut(float_series, num_bins, retbins=True, labels=False, duplicates='drop')

        logging.debug('Counting elements in each bin')
        bin_counts = pd.Series(float_series_binned).value_counts().sort_index()

        logging.info(f"Generated {len(bin_edges)-1} bins successfully")

        return bin_edges, bin_counts
    except Exception as e:
        logging.error(f"Error in interval_range: {e}", exc_info=True)
        return None, None

# main function
try:
    directory = "Document/Pickle/MutationRate/FloatList"
    documents = list_documents(directory)
    num_bins = 14
    counts = [0]*num_bins
    bin_edges = [0, 0.013, 0.02, 0.03, 0.041, 0.051, 0.062, 0.073, 0.083, 0.105, 0.117, 0.139, 0.186, 3.912]
    for chr in range(1,23):
        std_chr = f"chr{chr}"
        for input_fn in documents:
            if input_fn.split(sep="_")[1] == std_chr:
                pass
            else:
                continue
            logging.info(f"Progress started: {std_chr}")
            input_fn = os.path.join(directory, input_fn)
            with gzip.open(input_fn, 'rb') as f:
                adding_list = pickle.load(f)
            for rate in adding_list:
                for bin_num in range(num_bins):
                    if rate >= bin_edges[bin_num] and bin_num != 13:
                        continue
                    elif rate>= bin_edges[bin_num] and bin_num == 13:
                        counts[bin_num] +=1
                        break
                    else:
                        counts[bin_num-1] +=1
                        break
            logging.info("Progress is done")


    if bin_edges is not None and counts is not None:
        logging.info('Interval range determined')
        data_fn = f"roulette_mutation_rate_{num_bins}bins.txt"
        graph_fn = f"roulette_mutation_rate_sorted_{num_bins}bins.png"
        config_graph(num_bins)
        logging.info('Creating histogram')
        histogram(counts, num_bins, graph_fn)
        logging.info('Histogram created')
        logging.info('Saving data to txt file')
        data_save_txt(counts, bin_edges, data_fn)
        logging.info('Data saved to txt file')
    else:
        logging.error('Failed to determine interval range')
except Exception as e:
    logging.error(f'An error occurred: {e}', exc_info=True)