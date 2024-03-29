import pandas as pd
import pickle
import os
import numpy as np
import random
import sys
sys.path.append("E:/Alex_MM_tmp/salgomed_Alex/bin/Full_Scripts/sergii_code_to_process/")
sys.path.append("E:/Alex_MM_tmp/salgomed_Alex/bin/Full_Scripts/")
sys.path.append("E:/Alex_MM_tmp/salgomed_Alex/bin/")
sys.path.append("../scRegulation ")
import binomial_thm_enriched_list_nx as bn
import Correlation_Analysis as ca
from collections import Counter
import networkx as nx
import time
import dendro_code
from imp import reload
import plotting_08_02_2020_from_SD as pl_plot
reload(ca)
reload(dendro_code)
reload(pl_plot)
pcn = nx.read_edgelist("../data/PCN.txt")
#pcn = nx.read_edgelist("../../PCN.txt")

"""
Get Conservation, network enrichment and percent expression merged metrics 
for bootstrapped data.
If compare_dendro not given gets all metrics but conservation
but calculates dendrogram. 
If diff_expr != None data is assumed to be a dictionary of pre-computed correlation
where the keys are the sample names and the values are correltion data frames.
if need_bootstrap == False then diff_expr and data are already dictionaries
of bootstrapped data of differential expression.
Parameters:
data: Count file to find correlation and differential expression if diff_expr = None. 
      Dictionary of correlation data frames if diff_expr != None that are either
      bootstraps or samples to be bootstrapped if need_bootstrap = True.
      
save_dir: Directory to save output files to.
cell_list: List of endothelial cells if data is count file.
n_samps: Number of pseudo samples to generate if data is count data.
norm: Whether count data needs log scale normalization.
genes1: Genes to use in correlation, all other genes will be left out.
genes2: Target genes (receptors) to use in dendrogram and for final merged metric.
seed: Random seed to use.
itr: Number of bootstraps to genrate if data is not already bootstrapped.
compare_dendro: Dendrogram data frame where index is genes "Dendrogram" column is the 
                dendrogram position. If None analysis will stop after generating network
                enrichment, percent expression and the dendrogram for the current data.
                
markers: Genes of interest to label in plots in downstream analysis.
fraction: Minimum fraction of genes to be used in correlation.
diff_expr: Dictionary of precomputed differential expression data frames. Columns must include 'avg_logFC'
           to differentiate upregulated from downregulated genes. 'PCT.1' for using percent exprssion of 
           the genes as a metric in merge metric. 'Q-Val' (adjusted P-Value) for prioritizing genes for the network.
           
need_bootstrap: Whether or not the data needs bootstrapping or if its a dictionary of bootrstapped data frames.
Output:
gene_stats.pz: Pickle file of statistics for each bootstrap.
Corr_BS: Folder of correlation for each sample/pseudo sample.
Diff_BS: Folder of differential expression statistics for each sample/pseudo sample.
BS_Corr_Dict.pzz: Dictionary of the correlation matrix of all the bootstrapped samples.
Pseudo_Sample_Cells.csv: List of which endothelial cells are in each pseudo sample.
Peak_Count: Count of how many times each gene appears in the peak.
"""

def pseudo_peak_process(data,
                        save_dir,
                        cell_list = None,
                        n_samps = 10, 
                        genes1 = None,
                        genes2=None,
                        norm = True,
                        seed =1,
                        itr = 100,
                        compare_dendro=None,
                        markers = [],
                        fraction=.05, 
                        diff_expr = None,
                        need_bootstrap = True):
    #data = pd.read_csv(data_file, index_col = 0)
    
    
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
   
    
    if  diff_expr != None and not need_bootstrap:#isinstance(diff_expr, pd.DataFrame)
        gene_stats_dict = diff_expr
        bs_corr_dict = data
        
        
    else:
    
        if diff_expr != None: #isinstance(diff_expr,pd.DataFrame):
            corr_dict = data
            diff_dict = diff_expr
            corr_df = merge_sample_dict(corr_dict)
            diff_df = merge_sample_dict(diff_dict)
            print("Corr DF Shape", corr_df.shape)
            
            full_data_process(data_file=corr_df,
                                save_dir=save_dir,
                                cell_list=cell_list,
                                genes1 = genes1,
                                genes2=genes2,
                                norm = norm,
                                seed =seed,
                                compare_dendro=compare_dendro,
                                markers =markers,
                                fraction=fraction,
                                diff_expr =diff_df)
            bs = get_bootstrap(list(corr_dict.keys()),itr = itr, seed = 1)
            
        else:
            full_data_process(data_file=data,
                                save_dir=save_dir,
                                cell_list=cell_list,
                                genes1 = genes1,
                                genes2=genes2,
                                norm = norm,
                                seed =seed,
                                compare_dendro=compare_dendro,
                                markers =markers,
                                fraction=fraction)
            
            cols =list(cell_list)
            
            random.seed(seed)
            random.shuffle(cols)
        
            other_cells = set(data.columns) - set(cols)
            
            pseudo_samp_dict = {}
            pseudo_samp_df = pd.DataFrame(index = cols,columns = ["Samp"])
            #size = len(cols)/n_samps
        
            for i in range(n_samps):
                pseudo_samp_dict[i]=cols[i:len(cols):n_samps]
                pseudo_samp_df.loc[pseudo_samp_dict[i],"Samp"] = i
                bs = get_bootstrap(list(pseudo_samp_dict.keys()),itr = itr, seed = 1)
                
            pseudo_samp_df.to_csv(os.path.join(save_dir,"Pseudo_Sample_cells.csv"))
            
            
            corr_dir = os.path.join(save_dir,"Corr_BS")
            if not os.path.exists(corr_dir):
                os.mkdir(corr_dir)
            
            print("Computing Correlaiton")
            curr_time = time.time()
            corr_dict = {}
            for i in pseudo_samp_dict:
                corr_dict[i] = ca.get_corr(data[pseudo_samp_dict[i]],genes1,genes2, log_scale = norm)#,fraction=fraction)
                corr_dict[i].to_csv(os.path.join(corr_dir,"PS%i_Corr.csv"%i))
            print("Finished Correlation: %f"%(time.time()-curr_time))
        
            diff_dir = os.path.join(save_dir,"diff_BS")
            if not os.path.exists(diff_dir):
                os.mkdir(diff_dir)
            
            print("Computing Differential Expression")
            curr_time = time.time()
            diff_dict = {}
            for i in pseudo_samp_dict:
                diff_dict[i] = ca.get_diff(data[pseudo_samp_dict[i]], data[other_cells],log_scale = norm)
                diff_dict[i].to_csv(os.path.join(diff_dir,"PS%i_Diff.csv"%i))
            print("Finished Differential Expression: %f"%(time.time()-curr_time))
             
            bs = get_bootstrap(list(pseudo_samp_dict.keys()),itr = itr, seed = 1)
        
        
        
        print("Computing median Correlaiton")
        bs_corr_dict = {}
        curr_time = time.time()
        for i in range(len(bs)):
            ns = bs[i]
            curr_dict_list = [corr_dict[j] for j in  ns]
            bs_corr_dict[i] = pd.concat(curr_dict_list)
            bs_corr_dict[i] = bs_corr_dict[i].groupby(bs_corr_dict[i].index).median()
        print("Finished Median Correlation: %f"%(time.time()-curr_time))
        
   
        
        gene_stats_dict = {}
        print("Computing Median Stats")
        curr_time = time.time()
        for i in range(len(bs)):
            ns = bs[i]
            curr_dict_list = [diff_dict[j] for j in  ns]
            gene_stats_dict[i] = pd.concat(curr_dict_list)
            gene_stats_dict[i] = gene_stats_dict[i].groupby(gene_stats_dict[i].index).median()
            if markers != None:
                gene_stats_dict[i].loc[gene_stats_dict[i].index.isin(markers),"Marker"] =1
            else:
                gene_stats_dict[i]["Marker"] =0
        print("Finished Median Stats: %f"%(time.time()-curr_time))   
        
    dump(bs_corr_dict, os.path.join(save_dir,"BS_Corr_Dict.pz"))  
    print("Computing Binomial")
    curr_time = time.time()
    print(gene_stats_dict[0].columns)
    for k in gene_stats_dict:
        pval = gene_stats_dict[k]["Q-Val"]
        goi = pval.loc[pval<.01].sort_values().index 
        lfc = gene_stats_dict[k]['avg_logFC']
        lfc_goi = lfc.loc[lfc>0].index
        goi = set(goi).intersection(lfc_goi)
        goi = list(pval.loc[(pval<.01)&pval.index.isin(goi)].sort_values().head(1000).index)

        gene_stats_dict[k]=gene_stats_dict[k].join(bn.nx_binom(pcn,goi)[["Binomial_Prob"]], how = "left")
        gene_stats_dict[k]["Log_BN"] = -np.log10(gene_stats_dict[k]["Binomial_Prob"])
    print("Finished Binomial: %f"%(time.time()-curr_time))  
    
    print("Computing Dendrogram")
    curr_time = time.time()
    dendro_dict = get_dendro_dict(bs_corr_dict, markers)
    print("Finished Dendrogram: %f"%(time.time()-curr_time)) 
    
    for k in dendro_dict:
        gene_stats_dict[k] = gene_stats_dict[k].join(dendro_dict[k], how = "outer")
    
    if not isinstance(compare_dendro, pd.DataFrame) :
        dump(gene_stats_dict,os.path.join(save_dir,"gene_stats.pz"))
        return gene_stats_dict
        
    print("Computing Dendrogram Distance")
    curr_time = time.time()
    dendro_dist_dict = get_dendro_dist_dict(dendro_dict)
    print("Finished Dendrogram Distance: %f"%(time.time()-curr_time))
    compare_dendro = get_dendro_dist_dict({"Compare":compare_dendro})["Compare"]
    
    
    print(gene_stats_dict[0].columns)
    print("Computing Sliding Average")
    curr_time = time.time()
    coi_dict = {"Independent-3":["PCT.1","Log_BN","Conservation"],
                "All-4":["PCT.1","Log_BN","Conservation", "Marker"]}
    ws = 21
    for k in dendro_dist_dict:
        conservation = get_conservation(dendro_dist_dict[k],compare_dendro)
        gene_stats_dict[k] =  gene_stats_dict[k].join(conservation)
        for i in coi_dict:
            coi = coi_dict[i]
            
            curr_df = gene_stats_dict[k].copy()
            genes = curr_df["Dendrogram"].dropna().sort_values().index
            curr_win = merge_metrics(genes,[curr_df[coi]])
            curr_win = norm_max(sliding_window(genes,curr_win, window_size=ws))
           
            gene_stats_dict[k] = gene_stats_dict[k].join(curr_win, how = "outer")                                              
            gene_stats_dict[k].columns = list(gene_stats_dict[k].columns[:-1]) +[i+"_WS%i"%ws]
    dump(gene_stats_dict,os.path.join(save_dir,"gene_stats.pz"))
    
    print("Finished Sliding Average: %f"%(time.time()-curr_time))  
    print("Getting Peak")
    peak_count = []
    #print(gene_stats_dict[k].columns)
    
    gene_stats_dict,peak_count = get_window_peak(gene_stats_dict,ws = ws)
  
   
    peak_count.to_csv(os.path.join(save_dir,"Peak_Count.csv"))
    
    dump(gene_stats_dict,os.path.join(save_dir,"gene_stats.pz"))
    
    return  gene_stats_dict
    #genes = dendro_dict[k].sort_values("Dendrogram").index
 

"""
Get Conservation, network enrichment and percent expression merged metrics 
for the whole data set.
If compare_dendro not given gets all metrics but conservation
but calculates dendrogram. 
If diff_expr != None data is assumed to be a dictionary of pre-computed correlation
where the keys are the sample names and the values are correltion data frames.
if need_bootstrap == False then diff_expr and data are already dictionaries
Parameters:
      data: Count file to find correlation and differential expression if diff_expr = None. 
            Dictionary of correlation data frames if diff_expr != None that are either
            bootstraps or samples to be bootstrapped if need_bootstrap = True.
      save_dir: Directory to save output files to.
      cell_list: List of endothelial cells if data is count file.
      n_samps: Number of pseudo samples to generate if data is count data.
      norm: Whether count data needs log scale normalization.
      genes1: Genes to use in correlation, all other genes will be left out.
      genes2: Target genes (receptors) to use in dendrogram and for final merged metric.
      seed: Random seed to use.
      itr: Number of bootstraps to genrate if data is not already bootstrapped.
      compare_dendro: Dendrogram data frame where index is genes "Dendrogram" column is the 
                      dendrogram position. If None analysis will stop after generating network
                      enrichment, percent expression and the dendrogram for the current data.
      markers: Genes of interest to label in plots in downstream analysis.
      fraction: Minimum fraction of genes to be used in correlation.
      diff_expr: Dictionary of precomputed differential expression data frames. Columns must include 'avg_logFC'
                 to differentiate upregulated from downregulated genes. 'PCT.1' for using percent exprssion of 
                 the genes as a metric in merge metric. 'Q-Val' (adjusted P-Value) for prioritizing genes for the network.
Output:
gene_stats_Full: List of all statistics generated for each gene.
"""

def full_data_process(data_file,
                        save_dir,
                        cell_list,
                        genes1 = None,
                        genes2=None,
                        norm = True,
                        seed =1,
                        compare_dendro=None,
                        markers = [],
                        fraction = .05,
                        diff_expr = None,
                        debug = False):
    #data = pd.read_csv(data_file, index_col = 0)
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    
    data = data_file
    
    
    corr_dir = os.path.join(save_dir,"Corr_BS")
    if not os.path.exists(corr_dir):
        os.mkdir(corr_dir)
    
    if isinstance(diff_expr, pd.DataFrame):
        gene_stats_df = diff_expr
        corr_df = data
    else:
        cols =list(cell_list)
    
        random.seed(seed)
        random.shuffle(cols)
    
        other_cells = set(data.columns) - set(cols)
        print("Computing Correlaiton For Full")
        curr_time = time.time()
        corr_df = ca.get_corr(data[cell_list],genes1,genes2, log_scale = norm, fraction=fraction)
        
        print("Finished Correlation for Full: %f"%(time.time()-curr_time))

        if debug:
            print(corr_df.shape)

        diff_dir = os.path.join(save_dir,"diff_BS")
        if not os.path.exists(diff_dir):
            os.mkdir(diff_dir)

        if debug: print(((data[cell_list]!=0).astype(int).sum(axis=1)/float(data[cell_list].shape[1])).sort_values(ascending=False).head())

        print("Computing Differential Expression For Full")
        curr_time = time.time()
        gene_stats_df = ca.get_diff(data[cell_list], data[other_cells],log_scale = norm, debug = debug)
        if debug: print(gene_stats_df.max())
        gene_stats_df.to_csv(os.path.join(diff_dir,"Full_Diff.csv"))
        print("Finished Differential Expression: %f"%(time.time()-curr_time))
     
    corr_df.to_csv(os.path.join(corr_dir,"Full_Corr.csv"))
    
    print("Computing Median Stats")
    curr_time = time.time()
    if markers != None:
        gene_stats_df.loc[gene_stats_df.index.isin(markers),"Marker"] =1
    else:
        gene_stats_df["Marker"] =0
    print("Finished Median Stats: %f"%(time.time()-curr_time))   
        
    print("Computing Binomial")
    curr_time = time.time()
    print(gene_stats_df.columns)

    pval = gene_stats_df["Q-Val"]
    goi = pval.loc[pval<.01].sort_values().index 
    lfc = gene_stats_df['avg_logFC']
    lfc_goi = lfc.loc[lfc>0].index
    goi = set(goi).intersection(lfc_goi)
    goi = list(pval.loc[(pval<.01)&pval.index.isin(goi)].sort_values().head(1000).index)

    gene_stats_df =gene_stats_df.join(bn.nx_binom(pcn,goi)[["Binomial_Prob"]], how = "left")
    gene_stats_df["Log_BN"] = -np.log10(gene_stats_df["Binomial_Prob"])
    print("Finished Binomial: %f"%(time.time()-curr_time))  
    
    print("Computing Dendrogram")
    curr_time = time.time()
    dendro_df = get_dendro_dict({"Full":corr_df}, markers)["Full"]
    print("Finished Dendrogram: %f"%(time.time()-curr_time)) 
   
    gene_stats_df = gene_stats_df.join(dendro_df, how = "outer")
    
    if not isinstance(compare_dendro, pd.DataFrame) :
        gene_stats_df.to_csv(os.path.join(save_dir,"gene_stats_Full.csv"))
        return gene_stats_df
        
    print("Computing Dendrogram Distance")
    curr_time = time.time()
    dendro_dist_df = get_dendro_dist_dict({"Full":dendro_df})["Full"]
    compare_dendro = get_dendro_dist_dict({"Compare":compare_dendro})["Compare"]
    print("Finished Dendrogram Distance: %f"%(time.time()-curr_time))
    
    print(gene_stats_df.columns)
    print("Computing Sliding Average Full")
    curr_time = time.time()
    coi_dict = {"Independent-3":["PCT.1","Log_BN","Conservation"],
                "All-4":["PCT.1","Log_BN","Conservation", "Marker"]}

    conservation = get_conservation(dendro_dist_df,compare_dendro)
    gene_stats_df =  gene_stats_df.join(conservation)
    ws = 21
    for i in coi_dict:
        coi = coi_dict[i]
        
        curr_df = gene_stats_df.copy()
        genes = curr_df["Dendrogram"].dropna().sort_values().index
        curr_win = merge_metrics(genes,[curr_df[coi]])
        curr_win = norm_max(sliding_window(genes,curr_win, window_size=ws))
       
        gene_stats_df = gene_stats_df.join(curr_win, how = "outer")                                              
        gene_stats_df.columns = list(gene_stats_df.columns[:-1]) +[i+"_WS%i"%ws]
        
    print("Finished Sliding Average: %f"%(time.time()-curr_time))  
    
    gene_stats_df.to_csv(os.path.join(save_dir,"gene_stats_Full.csv"))
    
    return  gene_stats_df
  
def pre_process_batch(df_dict,
                      cell_dict,
                      genes1 = None,
                      genes2=None,
                      fraciont = .05,
                      norm = True):
    count = 0
    diff_dict = {}
    corr_dict = {}
    for k in df_dict:
        if count%10 == 9:
            print("%i/%i"%(count+1,len(df_dict)))
        count+=1
        
        data = df_dict[k]
        cell_list = cell_dict[k]
        other_cells = set(data.columns)-set(cell_list)
        
        curr_time = time.time()
        corr_df = ca.get_corr(data[cell_list],genes1,genes2, log_scale = norm, fraction=fraction)
        
        
        curr_time = time.time()
        gene_stats_df = ca.get_diff(data[cell_list], data[other_cells],log_scale = norm)
        
        df_dict[k] = gene_stats_df
        corr_dict[k] = corr_df
        
        
    return df_dict,corr_dict
        
        

    
"""
Finds genes located in a peak based on the first column of a data frame.
Parameters:
      df: Data frame to get peak from.
Output:
      Index of peak genes.
"""
def get_peak (df, thresh = .5):
    
    max_idx = df.idxmax()
    #if type(max_idx) == str:
    #    max_idx = [max_idx]
    max_idx = df.index.get_loc(max_idx[0])
    
    botidx = max_idx
    topidx = max_idx
    
    while botidx>0 and df.iloc[botidx,0] > .5:
        botidx -= 1
        
    while topidx<df.shape[0] and df.iloc[topidx,0] > .5:
        topidx += 1
        
    botidx += 1
    
    return df.iloc[botidx:topidx].index    


"""
Splits columns into random equally sized lists.
Pameteres:
      col_list: List of columns to split.
      
      n_parts: Number of parts to split columns into
      
      itr: Number of times to generate random partitions.
      
      seed: Random seed to set.
      
output:
      cross_valid_columns: List of all the partitions
"""
def get_cross_valid_cols(col_list,n_parts=10,itr = 1, seed = 1):
    cross_valid_columns=[]
    random.seed(seed)
    for i in range(itr):
        rand_cols = list(col_list)
        random.shuffle(rand_cols)
        cross_valid_columns += [list(set(rand_cols)-set(rand_cols[x::n_parts])) for x in range(n_parts)]
        
    return cross_valid_columns 
 
"""
Get bootstrap of columns list.
Parameter:
      col_list: Columns list to sample from.
      
      n:  Number of elements to select in each bootstrap. Default is length of col_list.
      
      itr: Number of bootstraps to perform.
      
      seed: Random seed to set
output:
      cross_valid_columns: List of bootstrapped Columns.
"""
def get_bootstrap(col_list,n=None,itr = 100, seed = 1):
    cross_valid_columns=[]
    random.seed(seed)
    if n == None:
        n = len(col_list)
    for i in range(itr):
        rand_cols = list(col_list)
        cross_valid_columns.append(np.random.choice(rand_cols,n))
        
    return cross_valid_columns 

"""
Compute data frame for each correlation data frame in med_corr_dict and return as a dictionary of dendrogram data frames.
Parameters:
      med_corr_dict: Dictionary of correlation data frames to compute dendrograms for.
      
      curr_goi: Genes to label in dendrogram plots produced.
      
      BS: Whether or not the data is bootstrapped or not.
"""
   
def get_dendro_dict(med_corr_dict, curr_goi, BS = False):
    #dendro_dict_cv = {}
    dendro_dict = {}
    for k in med_corr_dict:
    
        if BS:
            dendro_dict[k] = {}
        
        if not BS:
            dendro_dict[k] = dendro_code.dendro(med_corr_dict[k].fillna(0),                      
                                genesSubset=curr_goi,
                                inhLoc=18)
        
        else:
            for i in med_corr_dict_cv[k]:

                inloc = 18
                dendro_dict[k][i] = dendro_code.dendro(med_corr_dict_cv[k][i].fillna(0),                     
                                genesSubset=curr_goi,
                                inhLoc=18)
                plt.close()
    return dendro_dict  

"""
Get dendrogram distance for each dendrogram dataframe in the dictionary and return as a dictionary of 
distance dataframes.
Parameteres:
      
      dendro_dict: Dictionary of dendrograms to compute distances from.
      
      BS: Whether the dictionary is a dictionary of bootstrapped data or not.
      
Output:
      dendro_dist_dict: Dictionary of dendrogram distances.
"""
def get_dendro_dist_dict(dendro_dict, BS = False): 
    count = 0
    if not BS:
        dendro_dist_dict = {}
        for k in dendro_dict:
            print(k)
            curr_df = dendro_dict[k][["Dendrogram"]].dropna()
            dendro_dist_dict[k] = pd.DataFrame([np.abs(curr_df["Dendrogram"].values -x) for x in curr_df["Dendrogram"].values],
                                                index = curr_df.index,
                                                columns = curr_df.index)
            #for c1 in dendro_dict[k].index:
             #   for c2 in dendro_dict[k].index:
              #      dendro_dist_dict[k].loc[c1,c2] =  abs(dendro_dict[k].loc[c1,"Dendrogram"] - dendro_dict[k].loc[c2,"Dendrogram"]) 
                                                                  
                   
    else:
        for k in dendro_dict:
            if k != "Human":continue
            print(k)
            #if k != "Choroid": continue
            #dendro_dist_dict_cv[k] = pd.DataFrame()
            dendro_dist_dict[k] ={}
            count = 0
            for k2 in dendro_dict[k]:
                count +=1
                if count %10 ==9: print(count)
                dendro_dist_dict[k][k2] = []
                vals = dendro_dict[k][k2]["Dendrogram"].dropna().values
                genes = dendro_dict[k][k2]["Dendrogram"].dropna().index
                for i in range(len(vals)):
                    #if c1 != "KDR": continue
                    dendro_dist_dict[k][k2].append([])
                    for j in range(len(vals)):
                        dendro_dist_dict[k][k2][i].append(abs(vals[i] - vals[j]))
        
        dendro_dist_dict[k][k2] = pd.DataFrame(dendro_dist_dict[k][k2], index=genes,columns = genes)
        
    return dendro_dist_dict
"""
Get the conservation for each gene based on dendrogram distance.
Parameters:
      dd_df1: Dendrogram distance dataframe of interest.
      
      dd_df2: Dendrogram distance dataframe to compare it too.
      
      n: Number of closest genes to look at to determine conservation.
Output:
      df: Dataframe of conservation values.      
"""
def get_conservation(dd_df1, dd_df2,n =50):
    dd_con = pd.DataFrame()
    for c in dd_df1.columns:
        if c not in dd_df2.columns:
            dd_con.loc[c,"Conservation"] = None
        else:
            g1 = dd_df1.sort_values(c).head(n).index
            g2 = dd_df2.sort_values(c).head(n).index
            
            dd_con.loc[c,"Conservation"] = len(set(g1).intersection(g2))
    return dd_con
    
"""
Combine metrics by adding they're sum normalized values together into one.
parameteres:
      
      genes: Genes to use when merging metrics.
      
      metric_list: Data frame of metrics to merge.
      
output:
      
      df: Data frame of merged metric values.
"""
def merge_metrics (genes,metric_list):
    df = pd.DataFrame()
    for m in metric_list:
        m = m.reindex(genes)
        m = m/m.sum()
        
        df = df.join(m,how = "outer")
    df = df.mean(axis=1)
    return df
 
"""
Normalizes the columns so minimum is 0 and maximum is one.
Parameters:
      df: Data frame to normalize.
      
output:
      df: Normalized df.
"""
def norm_max(df):
    return (df-df.min())/(df.max()-df.min())
    
      
"""
Smooths data using a sliding window.
Parameters:
      gene_order: Order of genes in dendrogram to slide over.
      
      data: Data frame of data to average.
      
      window_size: size of sliding window.
      
      col: Column to use of data to use.
      
      na_fill: What to fill na values with.
      
      
Output:
      window_val: Values after window normalization.
"""
def sliding_window(gene_order, 
                   data, 
                   window_size = 21,
                   wrap = True,
                   col = None,
                   na_fill = 0):
    if col == None:
        col = "Windowed_Average"
    
    data = data.reindex(gene_order)
    data = data.fillna(na_fill)
    data = norm_max(data)#(data -data.min())/(data.max()-data.min())
    half = int(window_size/2)
    other_half = window_size - half
    
    data = list(data.values)
    window_val = []
    
    for i in range(len(data)):
        if i <half:
            if wrap:
                curr_data = data[:i+other_half] #+ data[i-half:] 
            else: 
                continue
        elif i+other_half > len(data):
            if wrap:
                curr_data = data[i-half:] #+ data[:i+other_half - len(data)]
            else:
                continue
        else:
            curr_data = data[i-half:i+other_half]
        
        window_val.append(np.nanmean(curr_data))
    window_val = pd.DataFrame(window_val, index = gene_order,columns = [col])
    
    return window_val
    
 

"""
Plot dendrogram and heatmap using data generated by full_data_process.
Parameters:
      data_dir: Directory where data is located.
      
      df_dir: Location tos save intermediate dataframe.
      
      plt_dir: Locataion to save figure.
      
      title: Title of the figure.
      
      markerBS_num: Bootstrap number to use.
      
      mark_loc: Location of marker in datafrmae.
"""
def plot_dedndro(data_dir,df_dir,plt_dir,title,markersBS_num =0, mark_loc = 0):
    
    gene_stats_full = pd.read_csv(os.path.join(data_dir,"gene_stats_Full.csv"),index_col = 0)
    gene_stats_full["BN"] = gene_stats_full["Log_BN"]
    gene_stats_full["Fraction"] = gene_stats_full["PCT.1"]
    gene_stats_full["DD-Conservation"] = gene_stats_full["Conservation"]
    markers = list(gene_stats_full.loc[gene_stats_full["Marker"]==1].index)
    
    bar_df = os.path.join(df_dir,title +".hdf")
    gene_stats_full.to_hdf(bar_df, key = "df")

    df = pd.read_csv(os.path.join(data_dir,"Corr_BS/Full_Corr.csv"), index_col=0)
   
    #M= df.values
    #fig = plt.figure(figsize=(8, 12))
    curr_plt = pl_plot.plotting(df,markers[mark_loc:],markers[:mark_loc],bar_df,plt_dir, suffix=title)

"""
Compute sliding window for different metrics.
Parameters:
      
      gene_stats_df: Dataframe with various gene stats to use.
      
      ws: window size to use for sliding window.
      
      coi_dict: Dictionary of columns to suse.
      
output:
      gene_stats_df: Data frame of gene stats with the merged metrics added.
"""
def get_window (gene_stats_df, 
                    ws = 21,
                    coi_dict = {"Independent-3":["PCT.1","Log_BN","Conservation"],
                    "All-4":["PCT.1","Log_BN","Conservation", "Marker"]}):


    for i in coi_dict:
        coi = coi_dict[i]
        
        curr_df = gene_stats_df.copy()
        genes = curr_df["Dendrogram"].dropna().sort_values().index
        curr_win = merge_metrics(genes,[curr_df[coi]])
        curr_win = norm_max(sliding_window(genes,curr_win, window_size=30))
       
        gene_stats_df = gene_stats_df.join(curr_win, how = "outer")                                              
        gene_stats_df.columns = list(gene_stats_df.columns[:-1]) +[i+"_WS%i"%ws]
    
    return gene_stats_df
"""
Get the merged statistics and adds it to each data frame in gene_stats_dict and returns the peak
count for each gene.
Parameters:
      gene_stats_dict: Dictionary of dataframes of gene statistics to use use for peak count.
      
      ws: Window size to use for sliding window smoothing.
      
      coi_dict: Dictionary where values are the columns to merge and keys are what to call the
                new merged data column.
"""
def get_window_peak(gene_stats_dict, 
                    ws = 21,
                    coi_dict = {"Independent-3":["PCT.1","Log_BN","Conservation"],
                    "All-4":["PCT.1","Log_BN","Conservation", "Marker"]},
                    compare_dendro = None):
   
    
    if isinstance(compare_dendro, pd.DataFrame):
        print("Computing all dendro_dist")
        dendro_dict = {k:gene_stats_dict[k][["Dendrogram"]].dropna() for k in gene_stats_dict}
        
        dendro_dist_dict = get_dendro_dist_dict(dendro_dict)
        #print("Finished Dendrogram Distance: %f"%(time.time()-curr_time))
        compare_dendro = get_dendro_dist_dict({"Compare":compare_dendro})["Compare"]
   
   
    for k in gene_stats_dict:
    
        if isinstance(compare_dendro, pd.DataFrame):
            conservation = get_conservation(dendro_dist_dict[k],compare_dendro)
            try:
                gene_stats_dict[k] =  gene_stats_dict[k].join(conservation)
            except:
                pass
        
        for i in coi_dict:
            coi = coi_dict[i]
            
            curr_df = gene_stats_dict[k].copy()
            genes = curr_df["Dendrogram"].dropna().sort_values().index
            curr_win = merge_metrics(genes,[curr_df[coi]])
            curr_win = norm_max(sliding_window(genes,curr_win, window_size=30))
           
            gene_stats_dict[k] = gene_stats_dict[k].join(curr_win, how = "outer")                                              
            gene_stats_dict[k].columns = list(gene_stats_dict[k].columns[:-1]) +[i+"_WS%i"%ws]
  
    

    peak_count = []
    
    peak_list = []
    for ck in coi_dict:
        ck = ck +"_WS%i"%ws
        peak_count =[]
        for k in gene_stats_dict:
            genes = gene_stats_dict[k]["Dendrogram"].dropna().sort_values().index
            peak_count += list(get_peak(gene_stats_dict[k].loc[genes,[ck]]))
            
        peak_count = Counter(peak_count)
        peak_count = pd.DataFrame(peak_count.values(),list(peak_count.keys()),[ck])
        peak_list.append(peak_count)
    
    peak_df = pd.concat(peak_list, axis=1)
    
        
    return gene_stats_dict, peak_df
   
def merge_sample_dict(sample_dict):
    
    merge_df = pd.concat(list(sample_dict.values()))
    merge_df = merge_df.groupby(merge_df.index).median()
    return merge_df
    
    
def read_sra (sra_dir,sra, endo_cell_df):
    endo_cell_list = []
    all_count = pd.DataFrame()
    for f in os.listdir(sra_dir):
        #print(f)
        if sra not in f: continue
        print(f)
        curr_df = pd.read_csv(os.path.join(sra_dir, f), index_col=0)
        
        srs = f.split("_")[1].split(".csv")[0]
        
        curr_endo = list(endo_cell_df[sra +"_" +srs].dropna())
        
        curr_endo = [x +"_" +srs if x in all_count.columns else x for x in curr_endo]
        
        endo_cell_list += curr_endo
        
        all_count = all_count.join(curr_df, how = "outer", rsuffix = srs)

    return all_count, endo_cell_list
    
      
def pre_process_batch(df_dict,
                      cell_dict,
                      genes1 = None,
                      genes2=None,
                      fraciont = .05,
                      norm = True):
    count = 0
    diff_dict = {}
    corr_dict = {}
    for k in df_dict:
        if count%10 == 9:
            print("%i/%i"%(count+1,len(df_dict)))
        count+=1
        
        data = df_dict[k]
        cell_list = cell_dict[k]
        other_cells = set(data.columns)-set(cell_list)
        
        curr_time = time.time()
        corr_df = ca.get_corr(data[cell_list],genes1,genes2, log_scale = norm, fraction=fraction)
        
        
        curr_time = time.time()
        gene_stats_df = ca.get_diff(data[cell_list], data[other_cells],log_scale = norm)
        
        df_dict[k] = gene_stats_df
        corr_dict[k] = corr_df
        
        
    return df_dict,corr_dict
        
    

