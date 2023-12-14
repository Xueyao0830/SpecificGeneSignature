#!/usr/bin/env python


# define a function to do the stats test
# this is the correct methods

# how to choose the alpha？

#Assert is running on python 3
import sys
assert sys.version_info.major>2#Assert is running on python 3

import pandas as pd
import numpy as np
import scanpy as sc
import csv
import warnings

import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import statsmodels.api as sm
from scipy.stats import (ttest_1samp, ttest_rel,ttest_ind,wilcoxon,t as t_dbn)
from statsmodels.stats.multitest import multipletests
import argparse,os


sc.settings.verbosity = 0 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
#Other parameters
sc.settings.set_figure_params(dpi=80)

#################################
#         Parse inputs          #                
#################################



parser = argparse.ArgumentParser() 
#input parameters    
parser.add_argument('-i', '--sample_input', default=None)
parser.add_argument('-c', '--csv_input', default=None)
parser.add_argument('-n', '--group_name_input', default=None)

args = parser.parse_args()

if args:
    #input parameters
    sample_input = args.sample_input
    csv_input = args.csv_input
    group_name_input = args.group_name_input




###########################
###########################
###########################
######   Methods   #######
###########################
###########################
###########################
def cal_wilcoxon_and_fdr(sample,gene_lists,group,alpha=0.05):
      
    warnings.filterwarnings("ignore", category=UserWarning, module="scipy.stats")
    
    csv_file_path = f'../analysis.d/result/csv/{sample_file_name}_{gene_list_file_name}_{group}_result.csv'
      
    df = {}
    
    var_list = sorted(sample.loc[group].unique().tolist())
    

    for gene_list in gene_lists:
      df[f"{gene_list}_pval"] = []
      df[f"{gene_list}_fdr"] = []
      df[f"{gene_list}_reject"] = []
      

      
      for var in var_list:
            pdx_n = sample.loc[:, sample.loc[group] == var]
            pdx_not = sample.loc[:, sample.loc[group] != var]
            
                    
            pdx_median_n = pdx_n.loc[gene_lists[gene_list]].median(axis=1)
            no_pdx_median_n = pdx_not.loc[gene_lists[gene_list]].median(axis=1)
            statistic, pvalue = wilcoxon(pdx_median_n, no_pdx_median_n, alternative='greater')
            df[f"{gene_list}_pval"].append(pvalue)

      
      reject, corrected_pvalues = multipletests(df[f"{gene_list}_pval"],alpha, method='fdr_bh')[0:2] # alpha is FWER
      df[f"{gene_list}_fdr"] = corrected_pvalues
      df[f"{gene_list}_reject"] = reject

    df = pd.DataFrame(df)

    # Set the index
    index = var_list
    df = df.set_index(pd.Index(index))

    # save the csv file

    df.to_csv(csv_file_path, index=True, header=True)
    
    return df
  
  
def plot_boxplot(sample, gene_lists,group):
      
      fig_file_path = f'../analysis.d/result/figure/{sample_file_name}_{gene_list_file_name}_{group}_result.pdf'

      
      num_gene_lists = len(gene_lists)
      
      fig,axes = plt.subplots(num_gene_lists,1,figsize=(6, 3*num_gene_lists))
      
      var_list = sorted(sample.loc[group].unique().tolist())
      
      df = {}
      
      for i, gene_list in enumerate(gene_lists):
            merge_median = []
            
            for var in var_list:
                  pdx_n = sample.loc[:,sample.loc[group] == var]
                  median_values = pdx_n.loc[gene_lists[gene_list]].median(axis=1)
                  df[f"{gene_list}_{var}_median"] = median_values
                  merge_median.append(median_values)
                  
            merge_median_df = pd.DataFrame(merge_median, index=var_list).T
            
            merge_median_df.boxplot(ax=axes[i])
            axes[i].set_xlabel('sample')
            axes[i].set_ylabel('log_expression')
            axes[i].set_title(f'The gene expression in {gene_list} list between different samples')
      plt.tight_layout()  # Adjust layout
      
      plt.savefig(fig_file_path)
      
      
      
##########################
##########################
##########################
######   Execute   #######
##########################
##########################
##########################
#sample
sample = pd.read_csv(sample_input,sep='\t',header=0)
sample = sample.set_index('#H:hugo')
sample = sample.iloc[:,1:]
#log(x+1)
sample_numeric = sample.iloc[9:,:].apply(pd.to_numeric, errors='coerce')
sample.iloc[9:,:] = np.log2(sample_numeric+1)

#gene_lists
gene_lists = pd.read_csv(str(csv_input))
gene_lists = gene_lists.to_dict(orient='list')
gene_lists = {key: [value for value in values if not pd.isna(value)] for key, values in gene_lists.items()}

#ignore if the gene is not in the sample
for gene_list in gene_lists:
    original_gene_list = gene_lists[gene_list]
    filtered_gene_list = [v for v in original_gene_list if v in sample.index]
    gene_lists[gene_list] = filtered_gene_list

#group
group = str(group_name_input)


#naming
sample_file_name = os.path.basename(sample_input)
gene_list_file_name = os.path.basename(csv_input)
sample_file_name = os.path.splitext(sample_file_name)[0]
gene_list_file_name = os.path.splitext(gene_list_file_name)[0]

cal_wilcoxon_and_fdr(sample,gene_lists,group)
plot_boxplot(sample, gene_lists,group)


####################