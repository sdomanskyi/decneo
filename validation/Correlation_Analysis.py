import sys
sys.path.append(".")
sys.path.append("../scRegulation/")
sys.path.append("/mnt/home/paterno1/anaconda3/lib/python3.7/site-packages/")
#sys.path.append("/mnt/home/paterno1/anaconda3/lib/python3.7/site-packages/tables")
#sys.path.append("../bin")
import pandas as pd
import os
import numpy as np
import PanglaoDBannotation as pldba
import sergii_io
import re
from scipy import stats
from scipy.spatial.distance import cdist
import random
import time


def get_all_corr (annot_file,
                  out_excel,
                  out_dir,
                  annot_dir,
                  count_dir,
                  gene_list1 = None,
                  gene_list2 = None, 
                  keyword= None,
                  column = "Cell type annotation",
                  index_column = False,
                  method = "pearson",
                  m2h_file = None,
                  log_scale = True,
                  ignore_error = False,
                  debug = False):
                  
    #Readd in annotation data
    annot_df = pd.read_excel(annot_file, index_col = [0,1,2])
    
    #get location of count files
    count_dir = os.path.join(count_dir,"var","www","html","SRA","SRA.final")
    
    
    #get annotations related to keywords
    if keyword != None:
        if type(keyword) == str:
            if index_column:
            
                annot_df = annot_df.loc[ annot_df.index.get_level_values(column) == keyword]
            else:
                annot_df = annot_df.loc[annot_df[column] == keyword]
        else:
            #keyword = re.compile("|".join(keyword)
            if index_column:
            
                annot_df = annot_df.loc[ annot_df.index.get_level_values(column).isin(keyword)]
            else:
                annot_df = annot_df.loc[annot_df[column].isin(keyword)]
    corr_dict = {}
    annot_df.to_csv("../Annot_run.csv")
    
    #read in mouse to human gene conversion
    if m2h_file != None:
        m2h = pd.read_csv(m2h_file)
        m2h = m2h.loc[~(m2h.index.isna() | m2h.index.duplicated(keep = "first"))]
    
    #get list of things already done
    processed = os.listdir(out_dir)
    
    for i in set(annot_df.index.get_level_values("SRS accession")):
        #extract SRA
        sra = pldba.getSRA(i,annot_df)
        #get clusters
        clust = set(annot_df.loc[annot_df.index.get_level_values("SRS accession")==i].index.get_level_values("Cluster index"))
        #get clusters and cells
        clust_df = pldba.extractPanglaoDBannotation(annot_dir, annot_df,sra,i,False)
        
        #skip if file already proecessd
        if sra + "_" + i +".csv" in processed:
            continue
        
        #get count file
        r_file = os.path.join(count_dir, '%s%s.sparse.RData.h5' % (sra, '_' + i if i!='notused' else ''))
        
        #read in count data
        try:
            count_df = sergii_io.readRDataFile(r_file)#,takeGeneSymbolOnly = True)
        except:
            
            sys.stderr.write("Failed to read %s"%r_file)
            if ignore_error:
                continue
            else:
                count_df = sergii_io.readRDataFile(r_file)
                #break
            
        #convert mouse to human
        if m2h_file != None and ("Mus musculus" in set(annot_df.loc[annot_df.index.get_level_values("SRS accession")==i,"Species"])):
            count_df.index = [ m2h.loc[m2h["MGI.symbol"] == x,"HGNC.symbol"].values[0] if x in set(m2h["MGI.symbol"]) else x for x in count_df.index]
            
        #get cells matching keyword
        coi = clust_df.loc[clust_df["cluster"].isin(clust)].index
        count_df = count_df[coi]
        
        #get correlation
        count_df = count_df.loc[~count_df.index.duplicated(keep = "first")]
        #curr_corr =  get_corr(count_df,gene_list1,gene_list2,method,log_scale)
        try:
            curr_corr =  get_corr(count_df,gene_list1,gene_list2,method,log_scale)
        except:
    
            sys.stderr.write("Failed to find correlation for %s"%r_file)
            if ignore_error:
                continue
            else:
                curr_corr =  get_corr(count_df,gene_list1,gene_list1,method,log_scale)
                break
    
        #curr_corr.to_csv("test_curr_corr.csv")
        
        
        name = sra + "_" + i
        curr_corr.to_csv(os.path.join(out_dir,name +".csv"))
        """
        try:
            for c in curr_corr.columns:
            
                if c in corr_dict:
                    corr_dict[c] = corr_dict[c].join(curr_corr[[c]], how = "outer")
                    corr_dict[c].columns = list(corr_dict[c])[:-1] + [name]
                    
                else:
                    corr_dict[c] = curr_corr[c]
                    corr_dict[c].columns = [name]
                    
                    
        except:
            
               
            sys.stderr.write("Failed to loop through correlation %s"%f)
            if ignore_error:
                continue
            else:
                break
        """
            
    #xl_file = pd.ExcelWriter(out_excel)
    
    #for k in corr_dict:
    #    corr_dict[k].to_excel(xl_file,sheet_name = k)
        
    #xl_file.close()
            
        
def get_all_diff_exp (annot_file,
                  out_excel,
                  out_dir,
                  annot_dir,
                  count_dir,
                  gene_list1 = None,
                  gene_list2 = None, 
                  keyword= None,
                  column = "Cell type annotation",
                  index_column = False,
                  method = "pearson",
                  m2h_file = None,
                  log_scale = True,
                  ignore_error = False,
                  debug = False,
                  reverse = True):
                  
    #Readd in annotation data
    annot_df = pd.read_excel(annot_file, index_col = [0,1,2])
    
    #get location of count files
    count_dir = os.path.join(count_dir,"var","www","html","SRA","SRA.final")
    
    
    #get annotations related to keywords
    if keyword != None:
        if type(keyword) == str:
            if index_column:
            
                annot_df = annot_df.loc[ annot_df.index.get_level_values(column) == keyword]
            else:
                annot_df = annot_df.loc[annot_df[column] == keyword]
        else:
            #keyword = re.compile("|".join(keyword)
            if index_column:
            
                annot_df = annot_df.loc[ annot_df.index.get_level_values(column).isin(keyword)]
            else:
                annot_df = annot_df.loc[annot_df[column].isin(keyword)]
    corr_dict = {}
    annot_df.to_csv("../Annot_run.csv")
    
    #read in mouse to human gene conversion
    if m2h_file != None:
        m2h = pd.read_csv(m2h_file)
        m2h = m2h.loc[~(m2h.index.isna() | m2h.index.duplicated(keep = "first"))]
    
    #get list of things already done
    processed = os.listdir(out_dir)
    vals = list(set(annot_df.index.get_level_values("SRS accession")))
    if reverse:
        vals.reverse()
    for i in vals:
        #extract SRA
        sra = pldba.getSRA(i,annot_df)
        #get clusters
        clust = set(annot_df.loc[annot_df.index.get_level_values("SRS accession")==i].index.get_level_values("Cluster index"))
        #get clusters and cells
        clust_df = pldba.extractPanglaoDBannotation(annot_dir, annot_df,sra,i,False)
        
        #skip if file already proecessd
        if sra + "_" + i +".csv" in processed:
            continue
        
        #get count file
        r_file = os.path.join(count_dir, '%s%s.sparse.RData.h5' % (sra, '_' + i if i!='notused' else ''))
        
        #read in count data
        try:
            count_df = sergii_io.readRDataFile(r_file)#,takeGeneSymbolOnly = True)
        except:
            
            sys.stderr.write("Failed to read %s"%r_file)
            if ignore_error:
                continue
            else:
                count_df = sergii_io.readRDataFile(r_file)
                #break
              
        #convert mouse to human
        if m2h_file != None and ("Mus musculus" in set(annot_df.loc[annot_df.index.get_level_values("SRS accession")==i,"Species"])):
            count_df.index = [ m2h.loc[m2h["MGI.symbol"] == x,"HGNC.symbol"].values[0] if x in set(m2h["MGI.symbol"]) else x for x in count_df.index]
            
        #get other cells
        count_df = count_df.loc[~count_df.index.duplicated(keep = "first")]
        coi = clust_df.loc[clust_df["cluster"].isin(clust)].index
        other_df = count_df.loc[:,~count_df.columns.isin(coi)]
        
        #get cells matching keyword
        count_df = count_df[coi]
        
     
        #get correlation
        
        #curr_corr =  get_corr(count_df,gene_list1,gene_list2,method,log_scale)
        try:
            curr_corr =  get_corr(count_df,other_df,log_scale)
        except:
    
            sys.stderr.write("Failed to find correlation for %s"%r_file)
            if ignore_error:
                continue
            else:
                curr_corr =  get_diff(count_df,other_df,log_scale)
                #break
    
        #curr_corr.to_csv("test_curr_corr.csv")
        
        
        name = sra + "_" + i
        curr_corr.to_csv(os.path.join(out_dir,name +".csv"))
        """
        try:
            for c in curr_corr.columns:
            
                if c in corr_dict:
                    corr_dict[c] = corr_dict[c].join(curr_corr[[c]], how = "outer")
                    corr_dict[c].columns = list(corr_dict[c])[:-1] + [name]
                    
                else:
                    corr_dict[c] = curr_corr[c]
                    corr_dict[c].columns = [name]
                    
                    
        except:
            
               
            sys.stderr.write("Failed to loop through correlation %s"%f)
            if ignore_error:
                continue
            else:
                break
        """
            
    #xl_file = pd.ExcelWriter(out_excel)
    
    #for k in corr_dict:
    #    corr_dict[k].to_excel(xl_file,sheet_name = k)
        
    #xl_file.close()
            
        
      
        
    

def get_corr (count_df, gene_list1 = None, gene_list2 = None, method = "pearson", log_scale = True,fraction = .05):

    if log_scale:
        count_df = get_log_scale(count_df)
    
    count_df = count_df.loc[((count_df!=0).astype(int).sum(axis =1)/float(count_df.shape[1])) >fraction]
    
    count_df = count_df.fillna(0)
    #count_df.to_csv("count_test.txt")
    if gene_list1 != None:
        count_df = count_df.loc[count_df.index.isin(gene_list1)]
        
    if gene_list2 != None:
        #count_small_df = count_df.loc[count_df.index.isin(gene_list2)].T
        gene_list2 = set(gene_list2).intersection(count_df.index)
        #corr = pd.DataFrame()
        corr = {}
        corr = 1 - cdist(count_df.loc[gene_list2].values,count_df.values,metric = "correlation")
        #for g2 in gene_list2:
            #corr[g2] = []
            #for g1 in count_df.index:
             #   corr[g2].append(stats.pearsonr(count_df.loc[g1].values,count_df.loc[g2].values)[0])
                
        #corr = pd.DataFrame(corr.values(),columns = (corr.keys()),index = count.index)
        corr = pd.DataFrame(corr,index = list(gene_list2),columns = count_df.index).T
        #count_small_df.to_csv("Test_samll_Count.csv")
        #count_df = count_df.T
        #corr = count_small_df.apply(lambda s: count_df.corrwith(s))
        #print(corr.head())
            
    if gene_list2 == None:
        corr = 1 - cdist(count_df.values,count_df.values,metric = "correlation")
        corr = pd.DataFrame(corr, count_df.index,count_df.index)
        #corr = count_df.T.corr(method = method)
        
    return corr
    
def get_diff (count_df,other_df, log_scale = True, debug= False, fillna=True):

    ct = time.time()
    if fillna:
        count_df = count_df.fillna(0)
        other_df = other_df.fillna(0)
    
    if log_scale:
        count_df = get_log_scale(count_df)
        other_df = get_log_scale(other_df)
    if debug: sys.stderr.write("Log Scale:%f"%(ct-time.time()))
    
    ct = time.time()
    ind = set(count_df.index).intersection(other_df.index)
    count_df = count_df.loc[ind]
    other_df = other_df.loc[ind]
    #other_df.index = count_df.index
    #print(count_df.head().index)
    #print(other_df.head().index)
    results_df = pd.DataFrame()
    results_df["avg_logFC"] = (count_df.loc[ind].mean(axis=1)).fillna(0)
    results_df["avg_logFC"] -=(other_df.loc[ind].mean(axis=1)).fillna(0)
    if debug: sys.stderr.write("Mean Diff:%f"%(ct-time.time()))
    if debug: print("Mean Diff:%f"%(ct-time.time()))
    
    ct = time.time()
    ttest_res = stats.ttest_ind(count_df.fillna(0).T.values,other_df.fillna(0).T.values)
    if debug: sys.stderr.write("T-Test:%f"%(ct-time.time()))
    if debug: print("T-Test:%f"%(ct-time.time()))
    #print(len(ttest_res), len(ttest_res[0]))
    ct = time.time()
    results_df["T-Val"] = ttest_res[0]
    results_df["P-Val"] = ttest_res[1]
    results_df["Q-Val"] = results_df["P-Val"]*results_df.shape[0]
    results_df["PCT.1"] = (count_df!=0).astype(int).sum(axis=1)/float(count_df.shape[1])
    results_df["PCT.2"] = (other_df!=0).astype(int).sum(axis=1)/float(other_df.shape[1])
    if debug: sys.stderr.write("Other Metrics:%f"%(ct-time.time()))
    if debug: print("Other Metrics:%f"%(ct-time.time()))
    return results_df

    
    
def get_log_scale(df):
    df = (df/df.sum(axis = 0))*10000
    df = np.log1p(df)
    return df

def merge_corr_df(dir, gene, merge_list = None):
    df_list = []
    for f in os.listdir(dir):
        if merge_list != None  and f not in merge_list: continue
        curr_df = pd.read_csv(os.path.join(dir,f), index_col = 0)
        if gene not in curr_df.columns: continue
        
        curr_df = curr_df[[gene]]
        curr_df.columns = [f.split(".")[0]]
        curr_df = curr_df.loc[~curr_df.index.duplicated(keep = "first")]
        df_list.append(curr_df)
       
    df =  pd.concat(df_list, join = "outer", axis = 1)
    return df
    
def merge_Many_corr_df(dir,outdir, gene_list, suffix = "", merge_list = None):
    df_dict = {}
    for f in os.listdir(dir):
        if merge_list != None  and f not in merge_list: continue
        curr_df = pd.read_csv(os.path.join(dir,f), index_col = 0)
        for gene in gene_list:
        
            if gene not in curr_df.columns: continue
            
            curr_gene_df = curr_df[[gene]]
            curr_gene_df.columns = [f.split(".")[0]]
            curr_gene_df = curr_gene_df.loc[~curr_gene_df.index.duplicated(keep = "first")]
            
            if gene not in df_dict:
                df_dict[gene] = []
            
            df_dict[gene].append(curr_gene_df)
            

    for k in df_dict:  
        df =  pd.concat(df_dict[k], join = "outer", axis = 1)
        df.to_csv(os.path.join(outdir,"%s_%s.csv"%(k,suffix)))
    #return df
    
def extract_cells (annot_file,
                  out_excel,
                  out_dir,
                  annot_dir,
                  count_dir,
                  gene_list1 = None,
                  gene_list2 = None, 
                  keyword= None,
                  column = "Cell type annotation",
                  index_column = False,
                  method = "pearson",
                  m2h_file = None,
                  log_scale = True,
                  ignore_error = False,
                  debug = False):
                  
    #Readd in annotation data
    annot_df = pd.read_excel(annot_file, index_col = [0,1,2])
    
    #get location of count files
    count_dir = os.path.join(count_dir,"var","www","html","SRA","SRA.final")
    
    
    #get annotations related to keywords
    if keyword != None:
        if type(keyword) == str:
            if index_column:
            
                annot_df = annot_df.loc[ annot_df.index.get_level_values(column) == keyword]
            else:
                annot_df = annot_df.loc[annot_df[column] == keyword]
        else:
            #keyword = re.compile("|".join(keyword)
            if index_column:
            
                annot_df = annot_df.loc[ annot_df.index.get_level_values(column).isin(keyword)]
            else:
                annot_df = annot_df.loc[annot_df[column].isin(keyword)]
    corr_dict = {}
    annot_df.to_csv("../Annot_run.csv")
    

    processed = os.listdir(out_dir)
    
    all_counts = {}
    
    for i in set(annot_df.index.get_level_values("SRS accession")):
        #extract SRA
        sra = pldba.getSRA(i,annot_df)
        #get clusters
        clust = set(annot_df.loc[annot_df.index.get_level_values("SRS accession")==i].index.get_level_values("Cluster index"))
        #get clusters and cells
        clust_df = pldba.extractPanglaoDBannotation(annot_dir, annot_df,sra,i,False)
        
        #skip if file already proecessd
        if sra + "_" + i +".csv" not in processed:
            continue
            
        coi = clust_df.loc[clust_df["cluster"].isin(clust)].index
        
        name = sra + "_" + i
        
        all_counts[name] = list(coi)
    
    
    all_counts = pd.DataFrame(all_counts.values(), all_counts.keys()).T
    
    all_counts.to_csv("../Endothelial_Cells.csv")
    
def get_all_corr_hdf (hdf_expr,
                  out_excel,
                  out_dir,
                  annot_dir,
                  count_dir,
                  gene_list1 = None,
                  gene_list2 = None, 
                  keyword= None,
                  column = "Cell type annotation",
                  index_column = False,
                  method = "pearson",
                  m2h_file = None,
                  log_scale = True,
                  ignore_error = False,
                  debug = False):
                  
    
    #get list of things already done
    processed = os.listdir(out_dir)
    
    for k in ["Mus musculus", "Homo sapiens"]:
    
        curr_spec_df = pd.read_hdf(hdf_expr, key = k, index_col = 0)
        
        for i in set( curr_spec_df.columns.get_level_values("batch")):
            
            #skip if file already proecessd
            if i+".csv" in processed:
                continue
            
       
            try:
                count_df = curr_spec_df.loc[:,curr_spec_df.columns.get_level_values("batch") ==i]
            except:
                
                sys.stderr.write("Failed to read %s"%r_file)
                if ignore_error:
                    continue
                else:
                    count_df = sergii_io.readRDataFile(r_file)
                    #break
                
          
            #get correlation
            count_df = count_df.loc[~count_df.index.duplicated(keep = "first")]
            #curr_corr =  get_corr(count_df,gene_list1,gene_list2,method,log_scale)
            try:
                curr_corr =  get_corr(count_df,gene_list1,gene_list2,method,log_scale)
            except:
        
                sys.stderr.write("Failed to find correlation for %s"%r_file)
                if ignore_error:
                    continue
                else:
                    curr_corr =  get_corr(count_df,gene_list1,gene_list1,method,log_scale)
                    break
        
            #curr_corr.to_csv("test_curr_corr.csv")
            
            
            name = i +"_" + k.split(" ")[0]
            curr_corr.to_csv(os.path.join(out_dir,name +".csv"))
     
def get_bootstrap(col_list,n=None,itr = 100, seed = 1):
    cross_valid_columns=[]
    random.seed(seed)
    if n == None:
        n = len(col_list)
    for i in range(itr):
        rand_cols = list(col_list)
        cross_valid_columns.append(np.random.choice(rand_cols,n))
        
    return cross_valid_columns 
    
    
def get_all_diff_exp_hdf (hdf_expr,
                  out_excel,
                  out_dir,
                  annot_dir,
                  count_dir,
                  gene_list1 = None,
                  gene_list2 = None, 
                  keyword= None,
                  column = "Cell type annotation",
                  index_column = False,
                  method = "pearson",
                  m2h_file = None,
                  log_scale = True,
                  ignore_error = False,
                  debug = True,
                  reverse = True
                  ):
    
    ct = time.time()
    #get location of count files
    count_dir = os.path.join(count_dir,"var","www","html","SRA","SRA.final")
    
    if debug: sys.stderr.write("Count Dir Setup:%f"%(ct - time.time()))
    if debug: print("Count Dir Setup:%f"%(ct - time.time()))
    ct = time.time()
    
    corr_dict = {}
    
    #read in mouse to human gene conversion
    if m2h_file != None:
        m2h = pd.read_csv(m2h_file)
        m2h = m2h.loc[~(m2h.index.isna() | m2h.index.duplicated(keep = "first"))]
    
    #get list of things already done
    processed = os.listdir(out_dir)
    
    if debug: sys.stderr.write("Other Setup:%f"%(ct-time.time()))
    if debug: print("Other Setup:%f"%(ct-time.time()))
    
    for k in ["Mus musculus", "Homo sapiens"]:
    
        ct = time.time()
        curr_spec_df = pd.read_hdf(hdf_expr, key = k, index_col = 0)
        if debug: sys.stderr.write("Loading HDF Count:%f"%(ct-time.time()))
        if debug: print("Loading HDF Count:%f"%(ct-time.time()))
        
        
        vals = list(set(curr_spec_df.columns.get_level_values("batch")))
        if reverse:
            vals.reverse()
            
        for i in vals:
            
            #skip if file already proecessd
            if i+".csv" in processed:
                continue
            
            
            try:
                count_df = curr_spec_df.loc[:,curr_spec_df.columns.get_level_values("batch") ==i]
            except:
                
                sys.stderr.write("Failed to read %s"%r_file)
                if ignore_error:
                    continue
                else:
                    count_df = curr_spec_df.loc[:,curr_spec_df.columns.get_level_values("batch") ==i]
                    count_df = sergii_io.readRDataFile(r_file)
                    break
    

            #extract SRA
            #sra = pldba.getSRA(i,annot_df)
            #get clusters
            #clust = set(annot_df.loc[annot_df.index.get_level_values("SRS accession")==i].index.get_level_values("Cluster index"))
            #get clusters and cells
            #clust_df = pldba.extractPanglaoDBannotation(annot_dir, annot_df,sra,i,False)
            
            #skip if file already proecessd
            if  i +".csv" in processed:
                continue
            
            
            ct = time.time()
            #get count file
            r_file = os.path.join(count_dir, '%s.sparse.RData.h5' % (i if i!='notused' else ''))
            
            #read in count data
            try:
                other_count_df = sergii_io.readRDataFile(r_file)#,takeGeneSymbolOnly = True)
            except:
                
                sys.stderr.write("Failed to read %s"%r_file)
                if ignore_error:
                    continue
                else:
                    other_count_df = sergii_io.readRDataFile(r_file)
                    #break
            if debug: sys.stderr.write("Read Other:%f"%(ct-time.time()))
            if debug: print("Read Other:%f"%(ct-time.time()))
            
            ct = time.time()
            #convert mouse to human
            if m2h_file != None and k == "Mus musculus":
                other_count_df.index = [ m2h.loc[m2h["MGI.symbol"] == x,"HGNC.symbol"].values[0] if x in set(m2h["MGI.symbol"]) else x for x in other_count_df.index]
                
            #get other cells
            other_count_df = other_count_df.loc[~other_count_df.index.duplicated(keep = "first")]
            #coi = clust_df.loc[clust_df["cluster"].isin(clust)].index
            other_count_df = other_count_df.loc[:,~other_count_df.columns.isin(count_df.columns)]
            
            #get cells matching keyword
            #count_df = count_df[coi]
            
         
            #get correlation
            
            #curr_corr =  get_corr(count_df,gene_list1,gene_list2,method,log_scale)
            other_count_df = get_log_scale(other_count_df)
            if debug: sys.stderr.write("Fix Other:%f"%(ct-time.time()))
            if debug: print("Fix Other:%f"%(ct-time.time()))
            
            ct = time.time()
            try:
                curr_corr =  get_diff(count_df,other_count_df,False, debug)#log_scale)
            except:
                sys.stderr.write("Failed to find correlation for %s"%r_file)
                if ignore_error:
                    continue
                else:
                    curr_corr =  get_diff(count_df,other_count_df,log_scale)
                    #break
            #if debug: print(curr_corr.head())
            
            if debug: sys.stderr.write("Get Diff:%f"%(ct-time.time()))
            if debug: print("Get Diff:%f"%(ct-time.time()))
            #curr_corr.to_csv("test_curr_corr.csv")
            
            ct = time.time()
            name =  i
            curr_corr.to_csv(os.path.join(out_dir,name +".csv"))
            
            if debug: print(os.path.join(out_dir,name +".csv"))
            
            if debug:sys.stderr.write("Saving:%f"%(ct-time.time()))
            if debug: print("Saving:%f"%(ct-time.time()))
