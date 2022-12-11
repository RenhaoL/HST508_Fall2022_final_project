# Hi-C class script
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore
import argparse
import itertools
from tqdm import tqdm
import random
import warnings

parser = argparse.ArgumentParser(description="TAD analysis")

# get inputs using parser
parser.add_argument('--data_path', type=str, default="../data/ENCODE_bulk_rna_seq.csv")
parser.add_argument('--gene_loc_path', type=str, default="../data/gene_locations.tsv")
parser.add_argument('--tad_path', type=str, default="../data/TAD_strong_boundary_start_end.csv")

args = parser.parse_args()

def extract_data(data_path,gene_loc_path,tad_path,hivar_pctl=None,excl_chrom=['chrM','chrX','chrY']):
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
    # filter gene by variance across samples/tissues if thus specified
    if type(hivar_pctl)==int or type(hivar_pctl)==float:
        if hivar_pctl>=0 and hivar_pctl<=100:
            sc_df_filtered = filter_genes_by_variance(sc_df_filtered,hivar_pctl)
    # get rid of chromosomes in exclusion list
    gloc_filtered = gloc[gloc['seqname'].isin(set(gloc.seqname).difference(excl_chrom))]
    tad_filtered = tad[tad['chrom'].isin(set(gloc.seqname).difference(excl_chrom))]
    # get intersecting genes in both
    gene_list = set(sc_df_filtered.index).intersection(gloc_filtered.index)
    sc_df_filtered = sc_df_filtered.loc[gene_list]
    gloc_filtered = gloc_filtered.loc[gene_list]
    return sc_df_filtered, gloc_filtered, tad_filtered

def filter_genes_by_variance(sc_df,percentile): # not tested
    """
    filters a TPM (genes x tissue !!!!) dataset by highest percentile of variance across tissue/samples

    IMPORTANT: do not feed in normalized data. That would make this pointless.
    """
    sc_df = sc_df.T
    gene_vars = sc_df.var()
    most_var_genes = (gene_vars >= np.percentile(gene_vars,percentile))
    return sc_df.loc[:,most_var_genes].T # inverts back to genes x samples/tissues

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

def slide_boundary(chr, start, end, num_iter=5, step_size=10000):
    """
    shift a TAD boundary left or right [num_iter] times while retaining the size of the TAD. 
    return a list of new boundaries and the corresponding distances from the original TAD.
    """
    chr_length = get_chr_lengths()[chr]
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

def high_low_corr_genes(corr_df, percentile=90, high=True):
    """
    if high, return all gene pairs greater than percentile
    if low, return all gene pairs less than percentile
    """
    corr_matrix = corr_df.to_numpy()
    values = corr_matrix[np.triu_indices_from(corr_matrix, 1)].flatten()
    threshold = np.percentile(values, percentile)
    all_pairs = get_gene_pairs(corr_df.index)
    correlated_gene_pairs = []
    for pair in all_pairs:
        gene1, gene2 = pair[0], pair[1]
        corr = corr_df.loc[gene1,gene2]
        if high:
            if corr > threshold and corr < 0.99:
                correlated_gene_pairs.append((gene1,gene2))
        else:
            if corr < threshold:
                correlated_gene_pairs.append((gene1,gene2))
    return correlated_gene_pairs

def calc_gene_dist(same_chrom_gene_pair,gene_loc):
    """
    Calculate gene distance assuming genes are on the same chromosome
    Midpoint distance
    """
    gene1 = same_chrom_gene_pair[0]
    gene2 = same_chrom_gene_pair[1]
    gene1_srt = gene_loc.loc[gene1]['start']
    gene1_end = gene_loc.loc[gene1]['end']
    gene2_srt = gene_loc.loc[gene2]['start']
    gene2_end = gene_loc.loc[gene2]['end']
    return abs((gene1_srt+gene1_end)/2-(gene2_srt+gene2_end)/2)

def gene_dist_correlation(same_chrom_gene_pair, gene_corr_matrix):
    """Extract the correlation between two genes on the same chromosome from the chomosome specific correlation matrix

    Returns:
        float: correlation between the two given genes
    """
    assert len(same_chrom_gene_pair) == 2, "Only accept two genes as input"
    gene1, gene2 = same_chrom_gene_pair[0], same_chrom_gene_pair[1]
    assert str(gene1) in list(gene_corr_matrix.columns) and str(gene2) in list(gene_corr_matrix.columns), "The input gene pairs does not locate on the same chromosome"
    return gene_corr_matrix.loc[gene1, gene2]
    
def genes_in_same_tad(gene_pair,tg_dict,return_false_if_same_genes=True):
    """
    given a pair of genes as a tuple, determine whether the genes are in the same TAD.
    """
    if return_false_if_same_genes: # whether to treat same gene as belonging to the same TAD
        if gene_pair[0]==gene_pair[1]:
            return False
    for v in tg_dict.values():
        if gene_pair[0] in v:
            if gene_pair[1] in v:
                return True
            else:
                return False
    return False

def corr_vs_dist(corr_df, gene_loc_df, percentile = 90, random_sample_size=1000, plot=False, title="Gene distance between high correlated and low correlated genes",save="../results/corr_vs_dist_plot.png", **kwargs):
    
    """compare the genetic distance vs. the gene correlation by random sample 1000 genes from all correlations.
        **kwargs will be passed down to pd.plot.hist()

    Returns:
        the high correlated and other genes pairs with correlation and genetic distance. 
    """
    
    # calculate the highly correlated genes
    high_corr_gene_pairs = high_low_corr_genes(corr_df, percentile=percentile)
    low_corr_gene_pairs = high_low_corr_genes(corr_df, percentile=percentile, high=False)
    
    high_corr_gene_pairs_df = pd.DataFrame(high_corr_gene_pairs, columns=["gene1", "gene2"])
    low_corr_gene_pairs_df  = pd.DataFrame(low_corr_gene_pairs, columns=["gene1", "gene2"])

    # random sample genes from each category (highly and other correlated)
    high_corr_gene_pairs_df_random = high_corr_gene_pairs_df.sample(n=random_sample_size, random_state=42)
    high_corr_gene_pairs_list_random = list(high_corr_gene_pairs_df_random.itertuples(index=False, name=None))

    low_corr_gene_pairs_df_random  = low_corr_gene_pairs_df.sample(n=random_sample_size, random_state=42)
    low_corr_gene_pairs_list_random = list(low_corr_gene_pairs_df_random.itertuples(index=False, name=None))

    # extract the distance and correlations
    high_corr_gene_dis_chr = []
    high_corr_gene_dis_corr = []

    for pair_of_genes in tqdm(high_corr_gene_pairs_list_random):
        high_corr_gene_dis_chr.append(calc_gene_dist(pair_of_genes, gene_loc_df))
        high_corr_gene_dis_corr.append(gene_dist_correlation(pair_of_genes, corr_df))
    
    low_corr_gene_dis_chr = []
    low_corr_gene_dis_corr = []

    for pair_of_genes in tqdm(low_corr_gene_pairs_list_random):
        low_corr_gene_dis_chr.append(calc_gene_dist(pair_of_genes, gene_loc_df))
        low_corr_gene_dis_corr.append(gene_dist_correlation(pair_of_genes, corr_df))
    
    # add the distance and correlation to the dataframe
    high_corr_gene_pairs_df_random["gene_dis"] = high_corr_gene_dis_chr
    high_corr_gene_pairs_df_random["gene_dis_corr"] = high_corr_gene_dis_corr
    high_corr_gene_pairs_df_random = high_corr_gene_pairs_df_random.reset_index().drop("index", axis=1)

    low_corr_gene_pairs_df_random["gene_dis"] = low_corr_gene_dis_chr
    low_corr_gene_pairs_df_random["gene_dis_corr"] = low_corr_gene_dis_corr
    low_corr_gene_pairs_df_random = low_corr_gene_pairs_df_random.reset_index().drop("index", axis=1)

    # for plotting
    if plot:
        combined_df = pd.concat([high_corr_gene_pairs_df_random[["gene_dis"]], low_corr_gene_pairs_df_random[["gene_dis"]]], axis=1)
        combined_df.columns = ["gene_dis_high_corr", "gene_dis_others"]
        combined_df.plot.hist(alpha=0.5, **kwargs)
        plt.title(title)
        plt.xlabel("genetic distance (bp)")
        plt.savefig(save)
        plt.show()
        
    return high_corr_gene_pairs_df_random, low_corr_gene_pairs_df_random

def get_gene_pairs(genes):
    """
    given a list of genes, return all possible pairs (not including self pairs)
    """
    return list(itertools.combinations(genes, 2))

def get_random_gene_pairs(genes, tg_dict, num_pairs=100):
    """
    given a list of genes, return a random list of gene pairs that are not located in TADs
    """
    gene_pairs = []
    count = 0
    while count < num_pairs:
        gene1 = random.choice(genes)
        gene2 = random.choice(genes)
        if not genes_in_same_tad((gene1, gene2), tg_dict) and gene1 != gene2:
            gene_pairs.append((gene1, gene2))
            count += 1
    return gene_pairs

def get_prop_hc_genes(gene_pairs, hc_gene_pairs):
    """
    return proportion of gene pairs that are highly correlated
    """
    hc_count = 0
    for pair in gene_pairs:
        if pair in hc_gene_pairs:
            hc_count += 1
    return hc_count / len(gene_pairs)

def hc_genes_in_tads(chromosome, all_genes_corr_df, gene_loc, tad, tg_dict, plot=True):
    """
    plot distribution of HC gene pair proportion in TAD vs random gene pairs
    """
    random.seed(10)
    hc_gene_pairs = high_low_corr_genes(all_genes_corr_df, 90)
    tad_prop_hc_genes = []
    random_prop_hc_genes = []

    for t in tad.index:
        if t not in tg_dict.keys():
            continue
        tad_chr, tad_start, tad_end = tad.loc[t,:]

        tad_gene_pairs = get_gene_pairs(get_genes_in_interval(tad_chr, tad_start, tad_end, gene_loc))
        if len(tad_gene_pairs) == 0:
            continue
        random_gene_pairs = get_random_gene_pairs(gene_loc.index, tg_dict, len(tad_gene_pairs))

        tad_prop_hc_genes.append(get_prop_hc_genes(tad_gene_pairs, hc_gene_pairs))
        random_prop_hc_genes.append(get_prop_hc_genes(random_gene_pairs, hc_gene_pairs))

    if plot:
        df = pd.DataFrame({"TAD": tad_prop_hc_genes, "random": random_prop_hc_genes})
        df.plot.hist(alpha=0.5)
        plt.xlabel("Proportion of HC gene pairs")
        plt.ylabel("Frequency")
        plt.title(chromosome)
        plt.legend(loc="upper right")
        plt.savefig(f"../results/hc_genes_in_tad_{chromosome}.png")
    
    return df
    
def sliding_windows(tpm, gene_loc, tad, tad_data, gene_loc_data, tad_size=30, plot=True):
    """
    plot lineplot and heatmap of sliding windows
    """
    tg_dict = tad_gene_dict(tad_data, gene_loc_data, tad_size)

    for t in tad.index:
        if t not in tg_dict.keys():
            continue
        tad_chr, tad_start, tad_end = tad.loc[t,:]
        tad_corr_df, tad_corr = calc_tad_coexp(tad_chr, tad_start, tad_end, gene_loc, tpm)
        new_boundaries, slide_distances = slide_boundary(tad_chr, tad_start, tad_end, 20)

        corr = [tad_corr]
        distances = [0] + slide_distances

        for new_chr, new_start, new_end in new_boundaries:
            new_corr_df, new_corr = calc_tad_coexp(new_chr, new_start, new_end, gene_loc, tpm)
            corr.append(new_corr)
        tad_location = "{}_{}-{}".format(tad_chr, int(tad_start), int(tad_end))

        if plot:
            # lineplot
            plt.plot(np.array(distances)/1000, corr)
            plt.title(tad_location)
            plt.xlabel("Distance from TAD boundary in kb")
            plt.ylabel("Average correlation")
            plt.savefig(f"../results/{tad_location}_sliding_window.png")
            # heatmap
            sns_plot = sns.clustermap(tad_corr_df)
            sns_plot.figure.savefig(f"../results/{tad_location}_heatmap.png")

def main():
    # ignore warnings
    warnings.filterwarnings("ignore")
    # read in data
    tpm_data, gene_loc_data, tad_data = extract_data(args.data_path, args.gene_loc_path, args.tad_path)
    # normalize tpm data
    norm_tpm = log2norm_tpm(tpm_data)
    # make dictionaries
    cg_dict = chromosome_gene_dict(gene_loc_data)
    tg_dict = tad_gene_dict(tad_data,gene_loc_data)

    chromosome_list = ['chr'+str(i+1) for i in range(19)]
    
    for chromosome in chromosome_list:
        print(f"Processing {chromosome}...")
        tpm, gene_loc, tad = get_genes_from_chromosome(chromosome,norm_tpm,tad_data,gene_loc_data)
        # analysis 1: correlation vs sliding windows
        sliding_windows(tpm, gene_loc, tad, tad_data, gene_loc_data, tad_size=30, plot=True)
        print("Analysis 1 done.")
        # analysis 2: are highly correlated genes in the same TAD?
        all_genes_corr_df = tpm.transpose().corr()
        analysis2_df = hc_genes_in_tads(chromosome, all_genes_corr_df, gene_loc, tad, tg_dict, plot=True)
        analysis2_df.to_csv(f"../results/analysis2_{chromosome}_df.csv")
        print("Analysis 2 done.")
        # analysis 3: correlation as a function of distance between genes
        analysis3_high_corr_df, analysis3_other_corr_df = corr_vs_dist(all_genes_corr_df, gene_loc_data, plot=True, title = f"Gene distance between high correlated and low correlated genes \n in {chromosome}", save=f"../results/corr_vs_dist_plot_{chromosome}.png")
        analysis3_high_corr_df.to_csv(f"../results/analysis3_{chromosome}_high_corr_df.csv")
        analysis3_other_corr_df.to_csv(f"../results/analysis3_{chromosome}_other_corr_df.csv")
        print(f"Analysis 3 done.")

if __name__ == "__main__":
    main()