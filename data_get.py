import gzip
import io
import numpy as np
import os, sys

def convert_to_float(line):
    converted_line = []
    for value in line:
        try:
            converted_line.append(float(value))
        except ValueError:
            converted_line.append(value)
    return converted_line

def load_data_set(file_path):
    formatted_data = []
    with io.BufferedReader(gzip.open(file_path, 'rb')) as data:
        for line in data:
            formatted_data.append(convert_to_float(line.decode('utf-8').strip().split('\t')))
    sample_names = formatted_data[0]
    cpg_sites = [line[0] for line in formatted_data[1:-1]]
    phenotypes = np.array(formatted_data[-1][1:])
    methylation_values = np.array([line[1:] for line in formatted_data[1:-1]])
    return sample_names, cpg_sites, phenotypes, methylation_values

def load_testdata(file_path):
    formatted_data = []
    with open(file_path,'r') as testfile:
        for line in testfile:
            formatted_data.append(convert_to_float(line.strip().split('\t')))
    sample_names = formatted_data[0]
    cpg_sites = [line[0] for line in formatted_data[2:]]
    phenotypes = np.array(formatted_data[1][1:])
    methylation_values = np.array([line[1:] for line in formatted_data[2:5609]])
    bad_samples = np.where(sum(np.isnan(methylation_values))>1000)[0]
    methylation_values = np.delete(methylation_values, bad_samples, axis = 1)
    methylation_values = methylation_values[~np.isnan(methylation_values).any(axis=1)]
    phenotypes = np.delete(phenotypes, bad_samples)
    return sample_names, cpg_sites, phenotypes, methylation_values

# def get_example_data():
#     test_data = load_data_set(f'./GSE74193_test.tsv.gz')
#     train_data = load_data_set(f'./GSE74193_train.tsv.gz')
#     return test_data, train_data

def get_testdata(file):
    test_data = load_testdata(file)
    return test_data
