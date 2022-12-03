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
parser.add_argument('--tad_path', type=str, default="../data/TAD_strong_boundary_start_end.csv")

args = parser.parse_args()

def extract_data(data_path,gene_loc_path,tad_path,excl_chrom=['chrM','chrX','chrY']):
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
    gloc = pd.read_csv(gene_loc_path,sep="\t",index_col=0)
    tad = pd.read_csv(tad_path)
    # get rid of genes with 0 exp across all samples
    sc_df_filtered = sc_df.loc[np.sum(sc_df,axis=1)!=0]
    # get rid of chromosomes in exclusion list
    gloc_filtered = gloc[gloc['seqname'].isin(set(gloc.seqname).difference(excl_chrom))]
    tad_filtered = tad[tad['chrom'].isin(set(gloc.seqname).difference(excl_chrom))]
    # get intersecting genes in both
    gene_list = set(sc_df_filtered.index).intersection(gloc_filtered.index)
    sc_df_filtered = sc_df_filtered.loc[gene_list]
    gloc_filtered = gloc_filtered.loc[gene_list]
    return sc_df_filtered, gloc_filtered, tad_filtered

def chromosome_gene_dict(gloc_data):
    """
    return dictionary where keys are chromosomes and values are gene lists
    """
    chromosomes = set(gloc_data.seqname)
    cg_dict = dict(zip(chromosomes, [None]*len(chromosomes)))
    for chrom in chromosomes:
        gloc = gloc_data[gloc_data.seqname==chrom]
        cg_dict[chrom] = list(gloc.index)
    return cg_dict

def get_genes_in_interval(chrom,start,end,gloc):
    """
    get genes in an interval on a chromosome
    """
    gloc_chr = gloc[gloc['seqname']==chrom]
    gloc_chr = gloc_chr[(gloc_chr['start'] >= start) & (gloc_chr['end'] < end)]
    return list(gloc_chr.index)

def tad_gene_dict(tad_locs,gloc,filter_bar = 3):
    tg_dict = dict(zip(list(range(len(tad_locs))),[[]]*len(tad_locs)))
    for i in range(len(tad_locs)):
        data = tad_locs.loc[i]
        tg_dict[i] = get_genes_in_interval(data['chrom'],data['start'],data['end'],gloc)
    if filter_bar > 0:
        tg_dict = {k:v for k,v in tg_dict.items() if len(v)>=filter_bar}
    return tg_dict

def log2norm_tpm(tpm_data):
    """
    returns log2-normalized tpm data
    """
    return zscore(np.log2(tpm_data+1),axis=1)

def get_genes_from_chromosome(chr_name,tpm_data,tad_data,gloc_data):
    """
    chr_name: e.g. 'chr1', the string corresponding to the chromosome you want to extract data on
    tpm_data: dataframe of tpms
    gene_loc_data: dataframe of gene locations
    ===
    returns filtered dataframes corresponding to chromosome of interest
    """
    gloc_filtered = gloc_data[gloc_data['seqname']==chr_name]
    tad_filtered = tad_data[tad_data['chrom']==chr_name]
    genes = gloc_filtered.index
    return tpm_data.loc[genes], gloc_filtered, tad_filtered

def get_chr_lengths(path_to_file="../data/chr_lengths"):
    """
    return a dict of chromosomes and their lengths
    """
    chr_lengths = {}
    infile = open(path_to_file)
    for line in infile:
        line = line.strip().split()
        chr = line[0]
        length = int(line[1])
        chr_lengths[chr] = length
    infile.close()
    return chr_lengths

def slide_boundary(chr, start, end, num_iter=5, x=0.2):
    """
    shift a TAD boundary left or right [num_iter] times with a step size of [x] * TAD size 
    while retaining the size of the TAD. return a list of new boundaries.
    """
    chr_length = get_chr_lengths()[chr]
    # set step size
    tad_length = end - start
    step_size = tad_length * x
    # TAD at start of chr -> shift right
    if start - (step_size * num_iter) < 0:
        direction = 1
    # TAD at end of chr -> shift left
    elif end + (step_size * num_iter) > chr_length:
        direction = 0
    # randomly shift left or right
    else:
        direction = np.random.randint(0,2)
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

def main():
    # read in data
    tpm_data, gene_loc_data, tad_data = extract_data(args.data_path, args.gene_loc_path, args.tad_path)
    # normalize tpm data
    norm_tpm = log2norm_tpm(tpm_data)
    # make dictionaries
    cg_dict = chromosome_gene_dict(gene_loc_data)
    tg_dict = tad_gene_dict(tad_data,gene_loc_data)

    chromosome_list = ['chr'+str(i+1) for i in range(19)]
    for chromosome in chromosome_list:
        tpm, gene_loc, tad = get_genes_from_chromosome(chromosome,norm_tpm,tad_data,gene_loc_data)
        # do stuff
        corr_df = tpm.transpose().corr() # correlation dataframe 
        avg_corr = corr_df.mean().mean()

if __name__ == "__main__":
    main()