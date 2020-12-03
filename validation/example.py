import sys
sys.path.append(".")
import Dendro_Dict as dd
import Correlation_Analysis as ca
import pandas as pd
import os
import PanglaoDBannotation as pldba
import sergii_io
import matplotlib
import time
import pickle
import h5py
matplotlib.use("Agg")

# Load in endothelial markers of interest and receptors
endo_goi = ['KDR','FLT1', 'FLT4', 'NRP1', 'NRP2', 'FGFR1', 'FGFR2','FGFR3', 'CXCR2', 'ROBO1','ROBO4', 'ENG',   
            'PDGFRA', 'PDGFRB', 'TEK', 'KIT', 'MET', 'CLEC14A', 'PNPLA2', 'CD36', 'CD47', 'VLDLR', 'PLXND1']

new_rec = "/mnt/home/paterno1/PL/receptorsListHugo_2555.txt"
all_recs =list(set( [x.strip() for x in open(new_rec).readlines()] + endo_goi))

## Run analysis on first species only 
out_dir = "results_hum/"
human_count_df = pd.read_csv("/mnt/home/paterno1/Michelle/Validation/human_count.csv", index_col = 0)

human_gene_stats_mk = dd.full_data_process(human_count_df,
                       out_dir,
                       None,
                       genes2=all_recs,
                       diff_expr = None,
                       markers = endo_goi)

## Run analysis on second species, allowing it to also calculate conservation 
out_dir = "results_mus/"
mouse_count_df = pd.read_csv("/mnt/home/paterno1/Michelle/Validation/mouse_count.csv", index_col = 0)
human_dendro = pd.read_csv("/mnt/home/paterno1/Michelle/Validation/results/gene_stats_Full.csv", index_col = 0)

mouse_gene_stats_mk = dd.full_data_process(mouse_count_df,
                       out_dir,
                       None,
                       genes2=all_recs,
                       diff_expr = None,
                       compare_dendro= human_dendro,
                       markers = endo_goi)