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

def slide_boundary(chr, start, end):
    """
    Given a TAD boundary, shift the boundary left or right while retaining the size of the TAD. 
    Return new boundary.
    """
    def get_new_boundary(start, end, step_size, direction): # 0 if left, 1 if right
        if direction == 1:
            new_start = start + step_size
            new_end = end + step_size
        else:
            new_start = start - step_size
            new_end = end - step_size
        return new_start, new_end
    # get length of chr
    infile = open("../data/chr_lengths")
    for line in infile:
        line = line.strip.split()
        if line[0] == chr:
            chr_length = line[1]
    infile.close()
    tad_length = end - start
    step_size = tad_length * 0.1
    # TAD at start of chr -> shift right
    if start - step_size < 0:
        new_start, new_end = get_new_boundary(start, end, step_size, 1)
    # TAD at end of chr -> shift left
    elif end + step_size > chr_length:
        new_start, new_end = get_new_boundary(start, end, step_size, 0)
    # randomly shift left (0) or right (1)
    else:
        new_start, new_end = get_new_boundary(start, end, step_size, np.random.randint(0,2))
    return new_start, new_end


def calc_tad_coexp():
    return None