# normalization.py
# HW2, Computational Genomics, Spring 2025
# andrewid: 

# WARNING: Do not change the file name; Autograder expects it.

import sys
import numpy as np
import matplotlib.pyplot as plt
from pprint import pprint
import os

PER_MILLION = 1/1000000
PER_KILOBASE = 1/1000


# Do not change this function signature
def rpkm(raw_counts, gene_lengths):
    """Find the normalized counts for raw_counts
    Returns: a matrix of same size as raw_counts
    """
    gene_count = len(raw_counts)
    samples_count = len(raw_counts[0])
    total_reads = np.sum(raw_counts, axis = 0)

    res = [[0] * samples_count for _ in range(gene_count)]
    for gene in range(len(raw_counts)):
        gene_length = gene_lengths[gene] #L(G)
        for sample in range(samples_count): 
            num_reads = raw_counts[gene][sample] #R(G)
            total_reads_samp = total_reads[sample] #N
            res[gene][sample] = (num_reads / gene_length) * ((10**9)/total_reads_samp)
    return res


   
# define any helper function here    

# Do not change this function signature
def size_factor(raw_counts):
    """Find the normalized counts for raw_counts
    Returns: a matrix of same size as raw_counts
    """
    num_genes = len(raw_counts)
    num_samples = len(raw_counts[0])
    res = [[0] * num_samples for _ in range(num_genes)]

    for sample in range(num_samples):
        k_vals = []
        for gene in range(num_genes):
            #calculate the ratios for each gene/sample pair
            observed_counts = raw_counts[gene][sample]
            product = np.prod(raw_counts[gene])
            geom_mean = product ** (1/num_samples)
            final_val = observed_counts / geom_mean
            k_vals.append(final_val)
        
        #get the size factor for each sample
        size_factor = np.median(k_vals) 

        #normalize the genes
        for gene in range(num_genes):
            observed_counts = raw_counts[gene][sample]
            res[gene][sample] = observed_counts / size_factor

    return res

# Assuming the functions `rpkm` and `size_factor` have been defined already

if __name__=="__main__":
    print("starting")
    raw_counts=np.loadtxt(sys.argv[1])
    gene_lengths=np.loadtxt(sys.argv[2])
    
    # print(raw_counts, gene_lengths)
    rpkm1=rpkm(raw_counts, gene_lengths) 

    size_factor1=size_factor(raw_counts)
    with open("size_factor_normalized_counts.txt", "w") as f:
        for i, row in enumerate(size_factor1):
            if i < len(size_factor1) - 1:
                f.write("\t".join(map(str, row)) + "\n")
            else:
                f.write("\t".join(map(str, row))) 
    print("done")
        
    # plot_boxplots(raw_counts, gene_lengths, output_dir='outputs')
    # print("done")

    # TODO: write plotting code here
    pass
