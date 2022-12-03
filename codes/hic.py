# Hi-C class script
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore

def extract_data(data_path="../data/single_cell_tpm.tsv",gene_loc_path="../data/gene_locations.tsv"):
    """
    return a dataframe with genes as rows and samples as columns
    make sure the genes are contained in the chromosome location info file (take intersection)
    """
    # read in data
    sc_df = pd.read_csv(data_path,sep="\t",index_col=0)
    # get rid of genes with 0 exp across all samples
    sc_df_filtered = sc_df.loc[np.sum(sc_df,axis=1)!=0]
    # read in gene location information
    gloc = pd.read_csv(gene_loc_path,sep="\t",index_col=0)
    return sc_df_filtered.loc[set(sc_df_filtered.index).intersection(gloc.index)]

def get_genes_from_chromosome(chromo_num):
    return None

def slide_boundary(chr, start, end, num_iter=5, x=0.2):
    """
    Given a TAD boundary, shift the boundary left or right [num_iter] times with a step size of [x] * TAD size 
    while retaining the size of the TAD. Return new boundaries.
    """
    def get_new_boundary(start, end, step_size, direction):
        new_boundaries = []
        for i in range(num_iter):
            if direction == 1: # shift right
                new_start, new_end = start + step_size, end + step_size
            else: # shift left
                new_start, new_end = start - step_size, end - step_size
            new_boundary = (chr, new_start, new_end)
            new_boundaries.append(new_boundary)
            start, end = new_start, new_end
        return new_boundaries
    # get length of chr
    infile = open("../data/chr_lengths")
    for line in infile:
        line = line.strip().split()
        if line[0] == chr:
            chr_length = int(line[1])
    infile.close()
    # set step size
    tad_length = end - start
    step_size = tad_length * x
    # TAD at start of chr -> shift right
    if start - (step_size * num_iter) < 0:
        new_boundaries = get_new_boundary(start, end, step_size, 1)
    # TAD at end of chr -> shift left
    elif end + (step_size * num_iter) > chr_length:
        new_boundaries = get_new_boundary(start, end, step_size, 0)
    # randomly shift left or right
    else:
        new_boundaries = get_new_boundary(start, end, step_size, np.random.randint(0,2))
    return new_boundaries


def calc_tad_coexp():
    return None