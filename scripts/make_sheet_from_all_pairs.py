import pandas as pd
import numpy as np
import sys

arg_list = sys.argv
lig = arg_list[1]
rec = arg_list[2]

end = '.h5.xlsx'
fname = 'v3_allpairs_0.0_choroid_'+lig+'_'+rec+end

data = pd.read_excel(fname)

f_comps = fname.split('_')
version = f_comps[0]
cutoff = f_comps[2]
ct1 = f_comps[4]
ct2 = f_comps[5].split('.')[0]
rl_pairs_unique = []

for nrow in range(data.shape[0]):
    row = data.iloc[nrow,]
    rl_pair_start = (row.ligand,row.receptor)
    #print(rl_pair_start)
    count_raw = row['count']
    rl_pairs_sub = 0
    row_pairs = row.pairs
  
    pairs = row_pairs.split(', ')
  
    for npair in range(len(pairs)):
        comps = pairs[npair].split('-')
        for comp in comps:
            if comp in rl_pair_start:
                rl_pairs_sub = rl_pairs_sub + 1
                break
            
    rl_pairs_unique.append(count_raw - rl_pairs_sub)

data['unique'] = rl_pairs_unique
data1 = data[['ligand','receptor','ram','count','count-1','unique','pairs']]

data1['perc_change'] = abs((data['unique']-data['count-1'])/data['count-1'])

data1.to_excel('intermediate_choroid/RL_pairs_unique_v3_choroid_'+cutoff+'_'+lig+'_'+rec+'.xlsx', index = False)
data1.to_csv('RL_pairs_unique_choroid_0.0_'+version+'_'+lig+'_'+rec+'.csv', index = False)