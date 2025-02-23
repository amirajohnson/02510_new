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
    final_k = 0
    significant_genes = []

    print(sorted_pvals)
    for i in range(len(pvals)):
        print(sorted_gene_names[i], sorted_pvals[i], (i+1)/n * alpha)
        if sorted_pvals[i] <= (i+1)/n * alpha:
            significant_genes.append(sorted_gene_names[i])

    print(significant_genes)
    return significant_genes



# define any helper function here    
# import matplotlib.pyplot as plt

#for part a
# def part_a():
#     # Load data
#     counts = np.loadtxt("size_factor_normalized_counts.txt", delimiter="\t")
#     labels = np.loadtxt("labels.txt", dtype=int)

#     print(len(counts))
#     # Split into case and control
#     treated = counts[:, labels == 1]
#     control = counts[:, labels == 2]

#     # Compute mean and dispersion for each gene
#     control_mean = np.mean(control, axis=1)
#     treated_mean = np.mean(treated, axis=1)

#     #compute std dev
#     control_std = np.std(control, axis=1)
#     treated_std = np.std(treated, axis=1)

#     #get disp
#     control_disp = control_std / (control_mean + 1e-10)
#     treated_disp = treated_std / (treated_mean + 1e-10)

#     #get log2 values
#     log2_control_mean = np.log2(control_mean + 1e-10)
#     log2_treated_mean = np.log2(treated_mean + 1e-10)

#     log2_control_disp = np.log2(control_disp + 1e-10)
#     log2_treated_disp = np.log2(treated_disp + 1e-10)

#     #for part b - log2 fold change
#     gene_names = np.loadtxt("GeneNames.txt", dtype=str)

#     log2_foldchange = np.log2((treated_mean) / control_mean)

#     #get the ten highest and ten lowest

#     lowest = np.argsort(log2_foldchange)[:10]
#     highest = np.argsort(log2_foldchange)[-10:]

#     genes = np.concatenate([lowest, highest])
#     assert(genes.shape == (20,)) #making sure we have ten highest and lowest

#     selected_gene_data = data[genes, :]

#     # Plotting the heatmap using matplotlib
#     fig, ax = plt.subplots(figsize=(10, 6))  # Adjust size as necessary

#     # Display the heatmap using imshow
#     cax = ax.imshow(selected_gene_data, cmap='coolwarm', aspect='auto')

#     # Set the tick labels for the y-axis (genes)
#     ax.set_yticks(np.arange(selected_gene_data.shape[0]))
#     ax.set_yticklabels(gene_names[genes])

#     # Set the tick labels for the x-axis (samples)
#     ax.set_xticks(np.arange(selected_gene_data.shape[1]))
#     # Here we assume that the number of columns in data matches the number of samples
#     # If gene_names is a list of sample names (not gene names), replace gene_names with it
#     ax.set_xticklabels([f"Sample {i+1}" for i in range(selected_gene_data.shape[1])], rotation=90, fontsize=10)


#     # Add color bar
#     fig.colorbar(cax, ax=ax, orientation='vertical', label='Expression Value')


#     # Add labels and title
#     plt.xlabel('Samples')
#     plt.ylabel('Genes (Top 10 Up & Down Regulated)')
#     plt.title('Gene Expression Heatmap for Top 10 Upregulated and Downregulated Genes')

#     plt.savefig("outputs/heatmap.png")

#     # Show the plot
#     plt.tight_layout()
#     plt.show()


    # # Plot scatter plot
    # plt.figure(figsize=(8, 6))
    # plt.scatter(log2_control_disp, log2_control_mean, color="red", alpha=0.5, label="Control")
    # plt.scatter(log2_treated_disp, log2_treated_mean, color="blue", alpha=0.5, label="Treated")
    # plt.xlabel("log2(Dispersion)")
    # plt.ylabel("log2(Mean)")
    # plt.title("Log2-Log2 Scatterplot of Gene Dispersion")
    # plt.legend()
    # # plt.grid(True)
    # plt.savefig("outputs/log2_log2_scatter_plot.png")
    # plt.show()

if __name__=="__main__":
    # Here is a free test case
    # genes=['a', 'b', 'c']
    # input1 = [0.01, 0.04, 0.1]
    # print(bh(genes, input1, 0.05))

    genes = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
    input = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    alpha = 0.1
    print(bh(genes, input, alpha))

    # data = np.loadtxt("size_factor_normalized_counts.txt", delimiter="\t", dtype = float)
    # gene_names = np.loadtxt("GeneNames.txt", dtype=str)  
    # labels_n = np.loadtxt("labels.txt", dtype=str)

    # print(data[0])
    # print(bh(gene_names, data[0], 0.05))

    #part_a() #for part a
    





