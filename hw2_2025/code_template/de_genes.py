# de_genes.py
# HW2, Computational Genomics, Spring 2025
# andrewid:

# WARNING: Do not change the file name; Autograder expects it.

import sys
import numpy as np
#from matplotlib import pyplot as plt

'''
Problem Description:
Determining whether the gene expressions in two conditions are statistically different
consists of rejecting the null hypothesis that the two data samples come from distributions with
equal means. To do this, we can calculate a p-value for each gene.
However. thresholding P-values to determine what fold changes are more significant than others
is not appropriate for this type of data analysis, due to the multiple hypothesis testing problem.
When performing a large number of simultaneous tests, the probability of getting a significant
result simply due to chance increases with the number of tests. In order to account for multiple
testing, perform a correction (or adjustment) of the P-values so that the probability of observing
at least one significant result due to chance remains below the desired significance level. We can
use The Benjamini-Hochberg (BH) adjustment.
The BH is defined as:
Let p1 ≤ ... ≤ pn be ordered p-values. Define
k = i : pi ≤ i
n α
and reject H1
0 ...Hk
0 . α is the false discovery rate(FDR) to be controlled. If no such i exists, then
no hypothesis will be rejected.
'''

# Do not change this function signature
def bh(gene_names, pvals, alpha):
    """(list, list, float) -> numpy array
    applies benjamini-hochberg procedure
    
    Parameters
    ----------
    genes: name of genes 
    pvalues: corresponding pvals
    alpha: desired false discovery rate
    
    Returns
    -------
    array containing gene names of significant genes.
    gene names do not need to be in any specific order.
    """
    #calculate the pvalue for each gene
    #perform a correction (or adjustment) of the P-values so that the probability of observing
    #at least one significant result due to chance remains below the desired significance level.
    #The BH is defined as:
    # Let p1 ≤ ... ≤ pn be ordered p-values. Define
    # k = i : pi ≤ i
    # n α
    # and reject H1
    # 0 ...Hk
    # 0 . α is the false discovery rate(FDR) to be controlled. If no such i exists, then
    # no hypothesis will be rejected.

    #get the ordered pvalues
    pvals = np.array(pvals)
    gene_names = np.array(gene_names)
    sorted_indices = np.argsort(pvals)
    sorted_pvals = pvals[sorted_indices]
    sorted_gene_names = gene_names[sorted_indices]
    n = len(pvals)
    final_k = -1
    significant_genes = []

    print(sorted_pvals)
    for i in range(len(pvals)):
        print(sorted_gene_names[i], sorted_pvals[i], (i+1)/n * alpha)
        if sorted_pvals[i] <= (i+1)/n * alpha:
            final_k = i
    
    if final_k != -1:
        significant_genes = sorted_gene_names[:final_k+1]
        return significant_genes.tolist()
    else: return []



# define any helper function here    

if __name__=="__main__":
    # Here is a free test case
    # genes=['a', 'b', 'c']
    # input1 = [0.01, 0.04, 0.1]
    # print(bh(genes, input1, 0.05))

    genes = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
    input = [0.001, 0.001, 0.01,  0.04, 0.1, 0.1, 0.1, 0.5]
    alpha = 0.001
    print(bh(genes, input, alpha))

    # data = np.loadtxt("size_factor_normalized_counts.txt", delimiter="\t", dtype = float)
    # gene_names = np.loadtxt("GeneNames.txt", dtype=str)  
    # labels_n = np.loadtxt("labels.txt", dtype=str)

    # print(data[0])
    # print(bh(gene_names, data[0], 0.05))

    #part_a() #for part a
    





