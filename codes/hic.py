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
parser.add_argument('--data_path', type=str, default="../data/ENCODE_bulk_rna_seq.csv")
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
    sc_df = pd.read_csv(data_path,index_col=0)
    sc_df = sc_df.loc[[idx for idx in sc_df.index if 'ENSMUSG' in idx]]
    gloc = pd.read_csv(gene_loc_path,sep="\t",index_col=0)
    tad = pd.read_csv(tad_path)
    # index manipulation
    sc_df_idx = [idx.split(".")[0] for idx in sc_df.index]
    sc_df.index = sc_df_idx
    gloc_idx = [idx.split(".")[0] for idx in gloc.index]
    gloc.index = gloc_idx
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

def tad_gene_dict(tad_locs,gloc,filter_bar=5):
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
    while retaining the size of the TAD. return a list of new boundaries and the
    corresponding distances from the original TAD.
    """
    chr_length = get_chr_lengths()[chr]
    # set step size
    tad_length = end - start
    step_size = 5000 # tad_length * x
    # default: shift right
    direction = 1
    # TAD at end: shift left
    if end + (step_size * num_iter) > chr_length:
        direction = -1
    new_boundaries = []
    distances = []
    step_size = direction * step_size
    distance_from_origin = step_size
    for i in range(num_iter):
        new_start, new_end = start + step_size, end + step_size
        new_boundaries.append((chr, new_start, new_end))
        distances.append(distance_from_origin)
        start, end = new_start, new_end
        distance_from_origin += step_size
    return new_boundaries, distances

def calc_tad_coexp(chr, start, end, gene_loc, tpm):
    """
    calculate the average pairwise correlation coefficient for a given genomic region
    """
    genes = get_genes_in_interval(chr, start, end, gene_loc)
    if len(genes) == 0:
        return None, 0
    tpm_subset = tpm.loc[genes,:]
    corr_df = tpm_subset.transpose().corr()
    avg_corr = np.mean(corr_df.to_numpy())
    return corr_df, avg_corr

def plot_corr_distance(distances, correlation):
    plt.plot(distances, correlation)
    plt.title("Average pairwise correlation")
    plt.xlabel("Distance from TAD boundary in kb")
    plt.ylabel("Avg correlation")
    plt.show()

def plot_corr_heatmap(corr_df):
    sns.heatmap(corr_df)
    plt.show()

# def get_highly_correlated_genes(corr_df, percentile=90):
#     """
#     given a correlation dataframe, return a list of the most highly correlated gene pairs.
#     """
#     values = corr_df.to_numpy().flatten()
#     threshold = np.percentile(values, 90)
#     genes = corr_df.index
#     pairs = []
#     for i in range(len(genes)):
#         gene1 = genes[i]
#         for j in range(len(genes)):
#             gene2 = genes[j]
#             corr = corr_df.iloc[i,j]
#             if corr > threshold:
#                 pairs.append((gene1, gene2))
#     return pairs


def genes_in_same_tad(gene_pair,tg_dict):
    """
    given a pair of genes as a tuple, determine whether the genes are in the same TAD.
    """
    gene1 = gene_pair[0]
    gene2 = gene_pair[1]
    for v in tg_dict.values():
        if gene1 in v:
            if gene2 in v:
                return True
            else:
                return False
    return False

def main():
    # read in data
    tpm_data, gene_loc_data, tad_data = extract_data(args.data_path, args.gene_loc_path, args.tad_path)
    # normalize tpm data
    norm_tpm = log2norm_tpm(tpm_data)
    # make dictionaries
    cg_dict = chromosome_gene_dict(gene_loc_data)
    tg_dict = tad_gene_dict(tad_data,gene_loc_data)

    chromosome_list = ['chr'+str(i+1) for i in range(19)]
    
    # analysis 1: correlation vs sliding windows
    for chromosome in chromosome_list:
        tpm, gene_loc, tad = get_genes_from_chromosome(chromosome,norm_tpm,tad_data,gene_loc_data)
        for t in tad.index:
            if t not in tg_dict.keys():
                continue
            tad_chr, tad_start, tad_end = tad.loc[t,:]
            tad_corr_df, tad_corr = calc_tad_coexp(tad_chr, tad_start, tad_end, gene_loc_data, tpm)
            # plot_corr_heatmap(tad_corr_df)
            new_boundaries, slide_distances = slide_boundary(tad_chr, tad_start, tad_end, 30)
            corr = [tad_corr]
            distances = [0] + slide_distances
            for new_chr, new_start, new_end in new_boundaries:
                new_corr_df, new_corr = calc_tad_coexp(new_chr, new_start, new_end, gene_loc_data, tpm)
                corr.append(new_corr)
            # plot_corr_distance(np.array(distances) / 1000 , corr)
    
    # # analysis 2: are highly correlated genes in the same TAD?
    # # get correlation matrix of all genes (genes * genes)
    # all_genes_corr_df = None
    # highly_correlated_genes = get_highly_correlated_genes(all_genes_corr_df)
    # pairs_in_tad = []
    # for pair in highly_correlated_genes:
    #     if genes_in_tad(pair):
    #         pairs_in_tad.append(pair)
    # # print frequency of pairs in the same TAD

    # analysis 3: correlation as a function of distance between genes

if __name__ == "__main__":
    main()