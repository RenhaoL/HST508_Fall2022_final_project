# Hi-C class script
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore
import argparse
import os

parser = argparse.ArgumentParser(description="TAD analysis")

# get inputs using parser
parser.add_argument('--data_path', type=str, default="../data/single_cell_tpm.tsv")
parser.add_argument('--gene_loc_path', type=str, default="../data/gene_locations.tsv")

args = parser.parse_args()

def extract_data(data_path,gene_loc_path,excl_chrom=['chrM','chrX','chrY']):
    """
    excl_chrom: list of strings corresponding to chromosomes to be excluded
    data_path: path to single cell tpm data
    gene_loc_path: path to file containing information on gene locations
    ===
    return a dataframe with genes as rows and samples as columns
    make sure the genes are contained in the chromosome location info file (take intersection)
    """
    # read in data
    sc_df = pd.read_csv(data_path,sep="\t",index_col=0)
    # get rid of genes with 0 exp across all samples
    sc_df_filtered = sc_df.loc[np.sum(sc_df,axis=1)!=0]
    # read in gene location information
    gloc = pd.read_csv(gene_loc_path,sep="\t",index_col=0)
    # get rid of gloc for chromosomes in exclusion list
    gloc_filtered = gloc[gloc['seqname'].isin(set(gloc.seqname).difference(excl_chrom))]
    # get intersecting genes in both
    gene_list = set(sc_df_filtered.index).intersection(gloc_filtered.index)
    return sc_df_filtered.loc[gene_list], gloc_filtered.loc[gene_list]

def log2norm_tpm(tpm_data):
    """
    returns log2-normalized tpm data
    """
    return zscore(np.log2(tpm_data+1),axis=1)

def get_genes_from_chromosome(chr_name,tpm_data,gloc_data):
    """
    chr_name: e.g. 'chr1', the string corresponding to the chromosome you want to extract data on
    tpm_data: dataframe of tpms
    gene_loc_data: dataframe of gene locations
    ===
    returns filtered dataframes corresponding to chromosome of interest
    """
    gloc_filtered = gloc_data[gloc_data['seqname']==chr_name]
    genes = gloc_filtered.index
    return tpm_data.loc[genes]

def slide_boundary():
    return None

def calc_tad_coexp():
    return None

def main():
    # read in data
    tpm_data, gene_loc_data = extract_data(args.data_path, args.gene_loc_path)
    
    # normalize tpm data
    norm_tpm = log2norm_tpm(tpm_data)

    chromosome_list = ['chr'+str(i+1) for i in range(19)]
    for chromosome in chromosome_list:
        tpm, gene_loc = get_genes_from_chromosome(chromosome,norm_tpm,gene_loc_data)
        # do stuff

if __name__ == "__main__":
    main()