#!/usr/bin/python

import sys
from scipy import stats
import pandas as pd
import networkx as nx



"""
Takes in a network x object and list of enriched genes and calculates the binomial
enrichment based on the number of enriched interaction there are.

Input:
    nx_obj: A networkx graph or edge list file.
    enr_gene: List of enriched genes.
    target_genes: List of target_genes. Default use all genes in background.
    background: List of genes to use as background probability. Default use all genes in the network.
    debug: Whether to print some debug statements
"""
def nx_binom(nx_obj,enr_gene,target_genes=False,background = False, debug = False):
    
    if not isinstance(nx_obj, nx.Graph):
        if isinstance(nx_obj, str):
            nx_obj = nx.read_edgelist(nx_obj)
            nx_obj = nx_obj.to_undirected()

        else:
            print("Error: nx_obj needs to be a networkx graph or a network edgelist file")
            return 0
    if background == False:
        background = nx_obj.nodes()
        
    if target_genes == False:
        target_genes = background
    
    
    binom_df = pd.DataFrame()

    enr_gene = set(enr_gene).intersection(background)
    
    e_e_interact = 0 #interact between 2 enr will be double counted - shoud divide by 2
    enr_interact = 0
    for e in enr_gene:
        eneigh = nx_obj.neighbors(e)
        enr_interact += len(set(eneigh) - set(enr_gene))
        e_e_interact += len(set(eneigh).intersection(enr_gene))
        
    #divide by 2 to not double count enriched interaction
    enr_interact += int(e_e_interact/2)
    
    #using networkx now so don't need to count edges (already read for me)
    p = float(enr_interact)/len(nx_obj.edges()) 
    
    if debug:
        print(enr_interact)
        print(len(nx_obj.edges()))
        print(p)
    
    trg_miss = 0
    trg_no_neigh = 0
    for g in target_genes:
        
        #keep track of target genes not in network for debugging
        if g not in nx_obj.nodes():
            trg_miss += 1
            continue
            
        #get gene list with background
        g_neigh = list(nx_obj.neighbors(g))
        g_hits = set(g_neigh).intersection(enr_gene)
        g_back = set(g_neigh).intersection(background)
        
        if len(g_back) == 0: 
            trg_no_neigh += 1
            continue
        
        #find binomial prob
        b_prob = stats.binom.sf(len(g_hits)-1,len(g_back), p)

        binom_df.loc[g, "Binomial_Prob"] = b_prob

        binom_df.loc[g, "Enriched_Interact"] = len(g_hits)
        binom_df.loc[g, "Total_Interact"] = len(g_back)
    
    if debug:
        print("%i Targets Missing"%trg_miss)
        print("%i Targts had no neighbors in background"%trg_no_neigh)
        
    return binom_df

"""
def nx_binom_carlo_method(nx_obj,enr_gene,target_genes=False,background = False, debug = False):
    
    if not isinstance(nx_obj, nx.Graph):
        if isinstance(nx_obj, str):
            nx_obj = nx.read_edgelist(nx_obj)
            nx_obj = nx_obj.to_undirected()

        else:
            print("Error: nx_obj needs to be a networkx graph or a network edgelist file")
            return 0
    if background == False:
        background = list(nx_obj.nodes())
        
    if target_genes == False:
        target_genes = background
    
    
    binom_df = pd.DataFrame()

    enr_gene = set(enr_gene).intersection(background)
    print("enr_gene", len(enr_gene))
    print(list(enr_gene)[:5])
    

    pval_list = []
    total_interact_list = []
    enrich_interact_list = []
    for g in list(nx_obj.nodes()):
        neigh = set(nx_obj.neighbors(g))
        total_interact = len(set(neigh).intersection(background))
        #enriched_interact = len(set(neigh).intersection(enr_gene).intersection(background))
        enriched_interact = len(set(neigh).intersection(enr_gene))
        
        #total_interact_list.append(total_interact)
        #enrich_interact_list.append(enriched_interact)
        
        #print(g, neigh,total_interact,enr_interact)
        
        
        if total_interact == 0:
            pval_list.append(0)
        else:
            pval_list.append(float(enriched_interact)/float(total_interact))
        
    #return(pval_list)
    #using networkx now so don't need to count edges (already read for me)
    p = sum(pval_list)/float(len(pval_list))

    if debug:
        print(len(nx_obj.edges()))
        print(p)
    
    trg_miss = 0
    trg_no_neigh = 0
    for g in target_genes:
    
        if g not in nx_obj.nodes():
            trg_miss += 1
            continue
        
        g_neigh = list(nx_obj.neighbors(g))
        g_hits = set(g_neigh).intersection(enr_gene)
        g_back = set(g_neigh).intersection(background)
        
        if len(g_back) == 0: 
            trg_no_neigh += 1
            continue
        
        b_prob = stats.binom.sf(len(g_hits)-1,len(g_back), p)

        binom_df.loc[g, "Binomial_Prob"] = b_prob

        binom_df.loc[g, "Enriched_Interact"] = len(g_hits)
        binom_df.loc[g, "Total_Interact"] = len(g_back)
    
    if debug:
        print("%i Targets Missing"%trg_miss)
        print("%i Targts had no neighbors in background"%trg_no_neigh)
        
    return binom_df
    
"""