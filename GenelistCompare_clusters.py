#Assert is running on python 3
import sys
assert sys.version_info.major>2#Assert is running on python 3

import pandas as pd
import numpy as np
import scanpy as sc
import csv

import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import statsmodels.api as sm
from scipy.stats import (ttest_1samp, ttest_rel,ttest_ind,wilcoxon,t as t_dbn)
from statsmodels.stats.multitest import multipletests
import os
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
parser.add_argument('-i', '--h5ad_input', default=None)
parser.add_argument('-c', '--csv_input', default=None)
parser.add_argument('-n', '--cluster_name_input', default=None)

args = parser.parse_args()

if args:
    #input parameters
    h5ad_input = args.h5ad_input
    csv_input = args.csv_input
    cluster_name_input = args.cluster_name_input



###########################
###########################
###########################
######   Methods   #######
###########################
###########################python
###########################

def cal_mannwhitneyu_fdr(sample,gene_lists,cluster_name='PAGODA_hc',alpha=0.05):
      '''
      input: 
      1. sample:h5ad file
      2. genelist_name: gene list, not a string
      3. genelist_score_name: gene list name, a string, format:''
      4. cluster_name: cluster name of annotation in h5ad file, a string. format:''
      
      '''
      
      
      csv_file_path = f'../analysis.d/result/csv/{sample_file_name}_{gene_list_file_name}_result.csv'
      
      
      # calculate the genelist scores
      
      gene_lists = gene_lists.to_dict(orient='list')
      gene_lists = {key: [value for value in values if not pd.isna(value)] for key, values in gene_lists.items()}
      
      
      ncluster = len(sample.obs[cluster_name].unique())
      df = {}
      annotation = sample.obs[cluster_name].unique().tolist()
      
      for gene_list in gene_lists:
            sc.tl.score_genes(sample,gene_lists[gene_list],score_name=str(gene_list),random_state=42)
            df[f"{gene_list}_pval"] = []  
            df[f"{gene_list}_FDR"] = []
            df[f"{gene_list}_reject"] =[]   

            for i in annotation:
                  cluster_i = sample.obs[sample.obs[cluster_name] == str(i)][str(gene_list)]
                  cluster_no_i = sample.obs[sample.obs[cluster_name] != str(i)][str(gene_list)]
                  statistic, pvalue = mannwhitneyu(cluster_i,cluster_no_i,alternative='greater',nan_policy='omit') 
                  
                  df[f"{gene_list}_pval"].append(pvalue)
      
            reject, pvals_corrected = multipletests(df[f"{gene_list}_pval"],alpha,method='hs')[0:2]

            df[f"{gene_list}_FDR"] = pvals_corrected
            df[f"{gene_list}_reject"] = reject

      df = pd.DataFrame(df)

      df.set_index(pd.Index(annotation), inplace=True)
      df.index.name = 'cluster'
      
      # save the csv file

      df.to_csv(csv_file_path, index=True)
      
      return df


##########################
##########################
##########################
######   Execute   #######
##########################
##########################
##########################

sample = sc.read(h5ad_input)
sample.var_names = [item.upper() for item in sample.var_names]
gene_lists = pd.read_csv(str(csv_input))
cluster_name = str(cluster_name_input)

sample_file_name = os.path.basename(h5ad_input)
gene_list_file_name = os.path.basename(csv_input)
sample_file_name = os.path.splitext(sample_file_name)[0]
gene_list_file_name = os.path.splitext(gene_list_file_name)[0]



cal_mannwhitneyu_fdr(sample,gene_lists,cluster_name)



####################
