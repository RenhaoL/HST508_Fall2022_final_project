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

def slide_boundary():
    return None

def calc_tad_coexp():
    return None