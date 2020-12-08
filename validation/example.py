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

# Get endothelial cell ids
endo_cell_df = pd.read_csv("/mnt/home/paterno1/Michelle/Validation/endo_cells.csv", index_col = 0)
endo_cells_list = list(endo_cell_df.columns)

## Run analysis on first species only, set hum_endo_list to list of endothelial cell ids in human_count_df
out_dir = "results_hum/"
human_count_df = pd.read_csv("/mnt/home/paterno1/Michelle/Validation/human_count.csv", index_col = 0)
hum_endo_list = [value for value in endo_cells_list if value in list(human_count_df.columns)] 

human_gene_stats_mk = dd.full_data_process(human_count_df,
                       out_dir,
                       hum_endo_list,
                       genes2=all_recs,
                       diff_expr = None,
                       markers = endo_goi)

## Run analysis on second species and set mus_endo_list to list of endothelial cell ids in mouse_count_df, allowing it to also calculate conservation 
out_dir = "results_mus/"
mouse_count_df = pd.read_csv("/mnt/home/paterno1/Michelle/Validation/mouse_count.csv", index_col = 0)
mus_endo_list = [value for value in endo_cells_list if value in list(mouse_count_df.columns)] 
human_dendro = pd.read_csv("/mnt/home/paterno1/Michelle/Validation/results_hum/gene_stats_Full.csv", index_col = 0)

mouse_gene_stats_mk = dd.full_data_process(mouse_count_df,
                       out_dir,
                       mus_endo_list,
                       genes2=all_recs,
                       diff_expr = None,
                       compare_dendro= human_dendro,
                       markers = endo_goi)