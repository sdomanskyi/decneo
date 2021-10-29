import matplotlib.pyplot as plt
import scipy
import scipy.stats
import collections
import os
import copy
import seaborn as sns
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans
import matplotlib.patheffects as path_effects
from matplotlib import colors as mcolors
def makePath(path):
    if not os.path.exists(path):
        os.makedirs(path)
        
    return True

def FindMax(ligs,recs):
    max_vals = []
    for ligand in ligs:
        for receptor in recs:
            uniquename = 'intermediate_choroid/choroid_heatmap_'+ligand+'_'+receptor+'_unique_0.0_values_nona.xlsx'
            unique = pd.read_excel(uniquename, index_col = 0)
    
            lig_rec_max = np.max(np.max(unique))
            max_vals.append(lig_rec_max)
    
    abs_max = np.max(max_vals)
    
    return abs_max


def heatmap(ligand,receptor):
    import openpyxl
    import pandas as pd
    import scipy
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    import os
    from matplotlib.patches import Rectangle




    sns.set(rc={'figure.figsize':(26.7,23.27)}) #figure Size
    sns.set(font_scale = 0.05)
    fname = 'RL_pairs_unique_v3_choroid_cutoff_0.0_'+ligand+'_'+receptor

    #'RL_pairs_unique_v3_choroid_cutoff_0.0_Macrophage_Endothelial'
    path = 'intermediate_choroid/'+fname + '.xlsx'
    #print(path)
    df = pd.read_excel(path,header = 0)
    #df1 = df.iloc[0:50,]
    #df2 = df1.pivot('ligand','receptor','count')

    param = fname.split('_')
    #species = param[6]
    ligand = param[7]
    receptor = param[8]
    cutoff = param[6]

    #C:\Users\bpham\Documents\Sanford\results 07 5 2021\cutoff_0.0\otherCellTypes_lig

    """
    sci_name = ''
    if species == 'human':
        sci_name = 'Homo sapiens. '
        elif species == 'mouse':
            sci_name = 'Mus musculus. ' 
    
    """

    l_path = 'results 09 17 2021 '+str(cutoff)+'/ligands/'+ ligand + '. dendrogram-heatmap-correlation-data.xlsx'
    r_path = 'results 09 17 2021 '+str(cutoff)+'/receptors/'+ receptor + '. dendrogram-heatmap-correlation-data.xlsx'




    l_dendo = pd.read_excel(l_path, header = 0)
    r_dendo = pd.read_excel(r_path,header = 0)


    ligand_order = l_dendo[l_dendo.columns[0]].unique()
    receptor_order = r_dendo[r_dendo.columns[0]].unique()


    gEC_choroidLigands = ['ADIPOQ', 'ANGPT1', 'ANGPT2', 'ANGPTL3', 'BMP2', 'BMP7', 'C3', 'C4B', 'DLL1',
                          'FN1', 'IL13', 'IL4', 'JAG1', 'JAG2', 'LIF', 'OSM', 'S1P', 'SEMA3A', 'SEMA3C',
                          'SEMA3E', 'SEMA4A', 'SLIT2', 'TGFB1', 'TGFB2', 'TGFB3', 'VEGFA'] # 26

    ligands_1777 = np.loadtxt('geneLists/ligands_1777.txt',dtype =str).tolist() #"genes of interest"
    igands_44 = np.loadtxt('geneLists/ligands_44.txt', dtype=str).tolist() #"known regulators"

    interesting_ligands =  gEC_choroidLigands # "known regulators" from geneslist in DECNEO # ligands_44 if panglao, gEC_choroidLigands if choroid
    interesting_ligands_in = []
    interesting_ligands_index = []
    for x in interesting_ligands:
        try:
            for z in range(0,len(ligand_order)):
                if x == ligand_order[z]:
                    interesting_ligands_index.append(z)
                    interesting_ligands_in.append(x)
        except:
            pass

    #print(interesting_ligands_index)
    #print(interesting_ligands_in)



    vals = ['unique']

    #print(ligand_order[interesting_ligands_index[0]])
    sig_ravgs_receptors = []
    sig_ravgs_ligands = []
    sig_ravgs_pairs = []
    row_pair = []
    col_pair = []
    for val in vals:
        df_final = df.pivot(index = 'ligand',columns = 'receptor',values = val).reindex(index = ligand_order, columns = receptor_order)
        #print(df_final.iloc[interesting_ligands_index[0],]).
        
        df_final.to_excel('intermediate_choroid/choroid_heatmap_'+ligand+'_'+receptor+'_'+val+'_'+cutoff+'_values.xlsx')
        df_final_nona = df_final.replace(np.nan, 0)
        df_final_nona.to_excel('intermediate_choroid/choroid_heatmap_' + ligand + '_' + receptor + '_' + val + '_' + cutoff + '_values_nona.xlsx')

   
    plt.figure()
    ligs = ['Fibroblast']
    recs = ['Endothelial']
    abs_max = FindMax(ligs,recs)
    ax = sns.heatmap(df_final,
                     cmap = sns.color_palette("rocket_r", as_cmap=True),
                     linewidths = 0.05,linecolor = None,
                     xticklabels = True, yticklabels = True,
                     square = True,vmin = 0)

    ax.xaxis.tick_top()  # x axis on top
    ax.xaxis.set_label_position('top')

    plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')

    for x in range(0,len(interesting_ligands_in)-1):
        an = ax.annotate(interesting_ligands_in[x],xy = (10,interesting_ligands_index[x] + 0.5),xycoords = 'data',textcoords = 'data',
                         xytext=(-150,interesting_ligands_index[x]+4), fontsize = 30,
                         arrowprops=dict(arrowstyle="fancy", color = 'black', alpha = 0.2))
        
    path = 'out_choroid/heatmaps/unique/'
    makePath(path)
    plt.savefig(path+'choroid_heatmap_'+ligand+'_'+receptor+'_'+val+'_'+cutoff+'.png',dpi = 1000)

    #plt.show()
  
    print('Done')
    return True


def heatmap_avg3sum(ligand,receptor, mode = 'sum'):
    import openpyxl
    import pandas as pd
    import scipy
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    import os
    from matplotlib.patches import Rectangle
    




    sns.set(rc={'figure.figsize':(26.7,23.27)}) #figure Size
    sns.set(font_scale = 0.05)
    fname = 'RL_pairs_unique_v3_choroid_cutoff_0.0_'+ligand+'_'+receptor


    path = 'intermediate_choroid/'+fname + '.xlsx'

    df = pd.read_excel(path,header = 0)


    param = fname.split('_')
    ligand = param[7]
    receptor = param[8]
    cutoff = param[6]


    l_path = 'results 09 17 2021 '+str(cutoff)+'/ligands/'+ ligand + '. dendrogram-heatmap-correlation-data.xlsx'
    r_path = 'results 09 17 2021 '+str(cutoff)+'/receptors/'+ receptor + '. dendrogram-heatmap-correlation-data.xlsx'



    l_dendo = pd.read_excel(l_path, header = 0,index_col = 0)
    r_dendo = pd.read_excel(r_path,header = 0,index_col = 0)
    
    
    ligand_order = l_dendo.index
    receptor_order = r_dendo.index
    
    lig_vals = l_dendo['Avg Combination of measures']
    lig_vals_max = np.max(lig_vals)
    lig_vals = lig_vals/lig_vals_max
    rec_vals = r_dendo['Avg Combination of measures']
    rec_vals_max = np.max(rec_vals)
    rec_vals = rec_vals/rec_vals_max
    lig_nrows = len(lig_vals)
    rec_nrows = len(rec_vals)
    
    array = np.zeros((lig_nrows,rec_nrows))
    
    for x in range(lig_nrows):
        for y in range(rec_nrows):
            if mode == 'sum':
                array[x,y] = (lig_vals[x] + rec_vals[y])/2 #Scale between 0 to 1
            elif mode == 'prod':
                array[x,y] = (lig_vals[x] * rec_vals[y])
            else:
                print('Not a valid mode')
                raise
    
    df= pd.DataFrame(array, index = ligand_order, columns = receptor_order)
    df.to_excel('intermediate_choroid/'+ligand+'_'+receptor+'_'+'avg3'+mode+'.xlsx')
    return True


import numpy as np
import pandas as pd
import scipy


def Unique_avg3sum_Corr(lig,rec,thresh,group=None,mode = 'sum', plot = False, corr = 'spearman'):
    import matplotlib.pyplot as plt
    
    parameters = {'axes.labelsize': 20,
          'axes.titlesize': 25}
    
    plt.rcParams.update(parameters)
    if group is not None:
        group = int(group)
        thresh = str(thresh)
    else:
        thresh = str(0.0)
    a_name = 'intermediate_choroid/choroid_heatmap_'+lig+'_'+rec+'_unique_0.0_values_nona.xlsx'
    if mode == 'prod':
        b_name = 'intermediate_choroid/'+lig+'_'+rec+'_avg3prod.xlsx'
    else:
        b_name = 'intermediate_choroid/'+lig+'_'+rec+'_avg3sum.xlsx'
    a = pd.read_excel(a_name,index_col = 0)
    b = pd.read_excel(b_name, index_col = 0)
    a_vals = []
    b_vals = []
    val_dict = {}
    if group is not None:
        annot_name = 'out_choroid/choroid'+'_'+lig+'_'+rec+'_annotation.xlsx'
        annot = pd.read_excel(annot_name, index_col = 0)
        annot = annot.loc[:,thresh+'_bounds'].dropna()
        interest = annot.iloc[group-1]
        lig_bounds = [int(a) for a in interest.split('_')[0].strip('][').split(', ')]
        rec_bounds = [int(a) for a in interest.split('_')[1].strip('][').split(', ')]
        for x in range(lig_bounds[0],lig_bounds[1] + 1):
            for y in range(rec_bounds[0],rec_bounds[1] + 1):
                a_vals.append(a.iloc[x,y])
                b_vals.append(b.iloc[x,y])
                
    else:
        bounds = a.shape
        for x in range(0,bounds[0]):
            for y in range(0,bounds[1]):
                a_vals.append(a.iloc[x,y])
                b_vals.append(b.iloc[x,y])
                
    for z in range(len(a_vals)):
        val_dict[a_vals[z]] = []
            
    for z in range(len(a_vals)):
        val_dict[a_vals[z]].append(b_vals[z])
            
    for z in range(len(a_vals)):
        val_dict[a_vals[z]] = np.mean(val_dict[a_vals[z]])
            
    val_dict = collections.OrderedDict(sorted(val_dict.items()))
    #rho, pval = scipy.stats.spearmanr(a_vals,b_vals,alternative = 'greater')
    if corr == 'spearman':
        rho, pval = scipy.stats.spearmanr(a_vals,b_vals)
    elif corr == 'pearson':
        rho, pval = scipy.stats.pearsonr(a_vals,b_vals)
    else:
        print('Invalid Correlation Mode')
        raise
    if plot == True:
        if group is None:
            g_label = 'All'
        else:
            g_label = str(group)
        
        
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
    #The p-value for a hypothesis test whose null hypotheisis
    #is that two sets of data are uncorrelated. 
    #See alternative above for alternative hypotheses.
    #pvalue has the same shape as correlation.
        ax1.set_title(mode+'_'+lig+'_'+rec+'_'+thresh+'_'+g_label+'_'+corr+' = ' + str(rho) + ' pval = ' +str(pval))
        ax1.set_xlabel('Number of Unique RL pairs' , fontsize = 25)
        ax1.set_xticks(np.arange(len(val_dict.keys())+1))
        ax1.set_yticks(np.arange(0,1.2,0.2))
        if  mode == 'sum':
            ax1.set_ylabel(mode + ' of % of max avg3combos / 2', fontsize = 25)
        else:
            ax1.set_ylabel(mode + ' of % of max avg3combos', fontsize = 25)
        ax1.plot(a_vals,b_vals,'ro')
        ax1.plot(val_dict.keys(),val_dict.values())
        
    
        ax1.xaxis.label.set_size(20)
        ax1.yaxis.label.set_size(20)
        
        plt.setp(ax1.get_xticklabels(), visible=True, ha='right', fontsize = 15)
        plt.setp(ax1.get_yticklabels(), visible=True, ha='right', fontsize = 15)
        fig.show()
        path = 'out_choroid/corrplots_avg3sum_vs_RLunique_spearman/'
        makePath(path)
        fig.savefig(path+lig+'_'+rec+'_'+thresh+'_'+g_label+'_'+corr+'_'+mode+'.png')
    return (rho,pval)




def ExportCorr(ligands,receptors,mode = 'sum'):
    if type(ligands) == str:
        ligands = [ligands]
    if type(receptors) == str:
        receptors = [receptors]
        
    df_dict = {}
    df_dict['Spearman'] = []
    df_dict['pval'] = []
    index = []
    for lig in ligands:
        for rec in receptors:
            index.append(lig+'_'+rec)
            out = Unique_avg3sum_Corr(lig,rec,0.2,group = None,mode = mode,plot = False)
            df_dict['Spearman'].append(out[0])
            df_dict['pval'].append(out[1])
    
            
            
    df = pd.DataFrame.from_dict(df_dict)
    df.index = index
    df.to_excel('intermediate_choroid/'+mode+'_avg3_RLunique_Correlation.xlsx')
    
    return df





def FindMaxAvg3(ligs,recs, mode = 'sum'):
    max_vals = []
    for ligand in ligs:
        for receptor in recs:
            
            avg3name = 'intermediate_choroid/'+ligand+'_'+receptor+'_'+'avg3'+mode+'.xlsx'
            avg3 = pd.read_excel(avg3name, index_col = 0)

            lig_rec_max = np.max(np.max(avg3))
            max_vals.append(lig_rec_max)
    
    abs_max = np.max(max_vals)
    
    return abs_max

def SumHeatmap(ligand,receptor, mode = 'sum'):
    import seaborn as sns
    
    sns.set(rc={'figure.figsize':(26.7,23.27)}) #figure Size
    sns.set(font_scale = 0.05)
    fname = 'intermediate_choroid/'+ligand+'_'+receptor+'_'+'avg3'+mode+'.xlsx'
    df = pd.read_excel(fname, index_col = 0)
    ligand_order = df.index.unique()
    receptor_order = df.columns.unique()
    #vals = ['count','count-1','unique','perc_change'] #'count','count-1',
    #df_interest = df[df['perc_change']>= 0.5]
    #interesting_ligands = df_interest['ligand'].unique()


    gEC_choroidLigands = ['ADIPOQ', 'ANGPT1', 'ANGPT2', 'ANGPTL3', 'BMP2', 'BMP7', 'C3', 'C4B', 'DLL1',
                          'FN1', 'IL13', 'IL4', 'JAG1', 'JAG2', 'LIF', 'OSM', 'S1P', 'SEMA3A', 'SEMA3C',
                          'SEMA3E', 'SEMA4A', 'SLIT2', 'TGFB1', 'TGFB2', 'TGFB3', 'VEGFA'] # 26

    ligands_1777 = np.loadtxt('geneLists/ligands_1777.txt',dtype =str).tolist() #"genes of interest"
    igands_44 = np.loadtxt('geneLists/ligands_44.txt', dtype=str).tolist() #"known regulators"

    interesting_ligands =  gEC_choroidLigands # "known regulators" from geneslist in DECNEO # ligands_44 if panglao, gEC_choroidLigands if choroid
    interesting_ligands_in = []
    interesting_ligands_index = []
    for x in interesting_ligands:
        try:
            for z in range(0,len(ligand_order)):
                if x == ligand_order[z]:
                    interesting_ligands_index.append(z)
                    interesting_ligands_in.append(x)
        except:
            pass
    ligs = ['Fibroblast','Macrophage','SMC','Pericyte','Endothelial']
    recs = ['Endothelial']
    abs_max = FindMaxAvg3(ligs,recs,mode = mode)
    plt.figure()
    ax = sns.heatmap(df,
                     cmap = sns.color_palette("rocket_r", as_cmap=True),
                     linewidths = 0.05,linecolor = None,
                     xticklabels = True, yticklabels = True,
                     square = True,vmin = 0, vmax = abs_max)

    ax.xaxis.tick_top()  # x axis on top
    ax.xaxis.set_label_position('top')

    plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')

    for x in range(0,len(interesting_ligands_in)-1):
        an = ax.annotate(interesting_ligands_in[x],xy = (10,interesting_ligands_index[x] + 0.5),xycoords = 'data',textcoords = 'data',
                         xytext=(-150,interesting_ligands_index[x]+4), fontsize = 30,
                         arrowprops=dict(arrowstyle="fancy", color = 'black', alpha = 0.2))
    path = 'out_choroid/heatmaps/avg3'+mode+'/'
    makePath(path)
    plt.savefig(path+'choroid_heatmap_'+ligand+'_'+receptor+'_avg3'+ mode +'.png',dpi = 1000)

    #plt.show()

    print('Done')
    return True



def Avg3andUnique_Heatmap(ligand,receptor, mode = 'sum', plot = True):
    sns.set(rc={'figure.figsize':(26.7,23.27)}) #figure Size
    sns.set(font_scale = 0.05)
    
    avg3name = 'intermediate_choroid/'+ligand+'_'+receptor+'_'+'avg3'+mode+'.xlsx'
    avg3 = pd.read_excel(avg3name, index_col = 0)
    
    uniquename = 'intermediate_choroid/choroid_heatmap_'+ligand+'_'+receptor+'_unique_0.0_values_nona.xlsx'
    unique = pd.read_excel(uniquename, index_col = 0)
    unique_max = np.max(np.max(unique))
    unique_scaled = unique/unique_max
    mean = (avg3+unique_scaled)/2 #scale from 0 to 1
    
    df = mean
    
    df.to_excel('intermediate_choroid/'+ligand+'_'+receptor+'_RLavg3mean_scaled.xlsx')

    if plot == True:
        ligand_order = df.index.unique()
        receptor_order = df.columns.unique()
        gEC_choroidLigands = ['ADIPOQ', 'ANGPT1', 'ANGPT2', 'ANGPTL3', 'BMP2', 'BMP7', 'C3', 'C4B', 'DLL1',
                              'FN1', 'IL13', 'IL4', 'JAG1', 'JAG2', 'LIF', 'OSM', 'S1P', 'SEMA3A', 'SEMA3C',
                              'SEMA3E', 'SEMA4A', 'SLIT2', 'TGFB1', 'TGFB2', 'TGFB3', 'VEGFA'] # 26

        ligands_1777 = np.loadtxt('geneLists/ligands_1777.txt',dtype =str).tolist() #"genes of interest"
        igands_44 = np.loadtxt('geneLists/ligands_44.txt', dtype=str).tolist() #"known regulators"

        interesting_ligands =  gEC_choroidLigands # "known regulators" from geneslist in DECNEO # ligands_44 if panglao, gEC_choroidLigands if choroid
        interesting_ligands_in = []
        interesting_ligands_index = []
        for x in interesting_ligands:
            try:
                for z in range(0,len(ligand_order)):
                    if x == ligand_order[z]:
                        interesting_ligands_index.append(z)
                        interesting_ligands_in.append(x)
            except:
                pass
        
    
        plt.figure()
        ax = sns.heatmap(df,
                         cmap = sns.color_palette("YlGnBu", as_cmap=True),
                         linewidths = 0.05,linecolor = None,
                         xticklabels = True, yticklabels = True,
                         square = True,vmin = 0,vmax = 1)

        ax.xaxis.tick_top()  # x axis on top
        ax.xaxis.set_label_position('top')

        plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')

        for x in range(0,len(interesting_ligands_in)-1):
            an = ax.annotate(interesting_ligands_in[x],xy = (10,interesting_ligands_index[x] + 0.5),xycoords = 'data',textcoords = 'data',
                             xytext=(-150,interesting_ligands_index[x]+4), fontsize = 30,
                             arrowprops=dict(arrowstyle="fancy", color = 'black', alpha = 0.2))
        path = 'out_choroid/heatmaps/RL_avg3'+ mode+'_mean/'
        makePath(path)
        plt.savefig(path+'choroid_heatmap_'+ligand+'_'+receptor+'_RLavg3'+ mode +'mean.png',dpi = 1000)

    
    return True

def MakeGroups(fname):
    file = pd.read_excel(fname, header = 0, index_col = 0)
    ligs = file.index.unique()
    recs = file.columns.unique()
    vals = []
    for nlig in range(len(ligs)):
        for nrecs in range(len(recs)):
            if file.iloc[nlig,nrecs] != 0:
                vals.append((nlig,nrecs))
    return (vals)

def groupNumbers2D(coors, minDist=10, minSize = 50):
    final_groups = []
    done = []
    
    tot = len(coors)
    counter = 0
    for coor in coors:
        temp = []
        x1 = coor[0]
        y1 = coor[1]
        counter += 1
        #print(str(counter) +' / '+str(tot))
        if coor in done:
            pass
        else:
            temp.append(coor)
            done.append(coor)
            for z in coors:
                x2 = z[0]
                y2 = z[1]
                diff_x = abs(x1-x2)
                diff_y = abs(y1-y2)
                if z != coor and z not in done and diff_x <= minDist and diff_y <= minDist:
                    temp.append(z)
                    done.append(z)
                    for member in temp:
                        for z in coors:
                            diff_x = abs(z[0]-member[0])
                            diff_y = abs(z[1] - member[1])
                            if z != member and z not in done and diff_x <= minDist and diff_y <= minDist:
                                temp.append(z)
                                done.append(z)
                    
                    final_groups.append(temp)
    y = []
    x = []
    groups = []
    counter = 0
    for group in final_groups:
        counter += 1
    
        for coor in group:
            y.append(coor[0])
            x.append(coor[1])
            groups.append(counter)
        
    df = pd.DataFrame({'lig' : y,
                  'rec' : x,
                  'group' : groups})
    
    
    return df

def MakeGroups2D(fname,minDist,minSizeClust):
    extend = fname.split('/')[1]
    param = extend.split('_')  #Macrophage_Endothelial_RLavg3mean_scaled
    lig = param[0]
    rec = param[1]
    a = MakeGroups(fname)
    b = groupNumbers2D(a, minDist = minDist, minSize = minSizeClust)
    b.to_hdf('2D_Dict.h5',key = lig+'_'+rec+'_minDist_'+str(minDist)+'_minSizeClust_'+str(minSizeClust), mode = 'a')
    return b



def MakeBounds(df,lig,rec, minSize = 50,  getTop = 3, findMax = False):
    fname = 'intermediate_choroid/'+lig+'_'+rec+'_RLavg3mean_scaled_nona.xlsx'
    data = pd.read_excel(fname,index_col = 0, header = 0)
    groups = df['group'].unique().tolist()
    groups_interest = []
    for group in groups:
        n_entries = df[df['group'] == group].shape[0]
        if n_entries >= minSize:
            groups_interest.append(group)
    bounds = dict()
    lig_mins = []
    lig_maxs = []
    rec_mins = []
    rec_maxs = []
    areas = []
    cluster_n = []
    counter = 1
    max_vals = []
    for group in groups_interest:
        df_interest = df[df['group'] == group]
        lig_min = np.min(df_interest['lig'].values)
        lig_max = np.max(df_interest['lig'].values)
        
        lig_mins.append(lig_min)
        lig_maxs.append(lig_max)

        rec_min = np.min(df_interest['rec'].values)
        rec_max = np.max(df_interest['rec'].values)
        
        rec_mins.append(rec_min)
        rec_maxs.append(rec_max)
        
        length = len(range(rec_min,rec_max+1))
        width = len(range(lig_min,lig_max+1))
        area = length * width
        areas.append(area)
        data_interest = data.iloc[lig_min:lig_max+1,rec_min:rec_max+1]
        data_interest_maxval = np.max(np.max(data_interest))
        max_vals.append(data_interest_maxval)
        
    for z in range(len(areas)):
        cluster_n.append(counter)
        counter += 1

    df_bounds = pd.DataFrame({'rec_min' : rec_mins,
                            'rec_max' : rec_maxs,
                            'lig_min' : lig_mins,
                            'lig_max' : lig_maxs,
                             'area' : areas,
                              'max_val' : max_vals}, index = groups_interest)
    df_bounds = df_bounds.sort_values(by='max_val',ascending = False)
    df_bounds.index = cluster_n
    if getTop != None:
        df_bounds = df_bounds.iloc[0:getTop,]
    if findMax == True:
        df_bounds = df_bounds.iloc[0:1,]
    return df_bounds


def GetRecsandLigs(fname):
    file = pd.read_excel(fname, index_col = 0, header = 0)
    ligs = file.index.unique()
    recs = file.columns.unique()
    return (recs,ligs)


def GetGroups(df,fname,out):
    recs_ligs = GetRecsandLigs(fname)
    params = fname.split('/')[1].split('_')
    lig_name = params[0]
    rec_name = params[1]
    nrows = df.shape[0]
    groups_ligs = dict()
    groups_recs = dict()
    nligs = []
    nrecs = []
    for z in range(nrows):
        row = df.iloc[z,]
        recs = []
        ligs = []
        for a in range(int(row.rec_min),int(row.rec_max)+1):
            recs.append(recs_ligs[0][a])
        for b in range(int(row.lig_min),int(row.lig_max)+1):
            ligs.append(recs_ligs[1][b])
        nligs.append(len(ligs))
        nrecs.append(len(recs))
        ligs = str(ligs).strip('][').replace("'",'')
        recs = str(recs).strip('][').replace("'",'')
        
        groups_ligs[df.index[z]] = ligs
        groups_recs[df.index[z]] = recs
        
    ligs = pd.Series(groups_ligs)
    recs = pd.Series(groups_recs)
    nligs = pd.Series(nligs, index = df.index)
    nrecs = pd.Series(nrecs, index = df.index)
    df_annot = pd.DataFrame({'recs' : recs, 'ligs' : ligs, 'n_recs' : nrecs, 'n_ligs' : nligs})
    fout = out+'Top_borders/'
    makePath(fout)
    df_annot.to_excel(fout +lig_name+'_'+rec_name+'_top_borders.xlsx')
    return df_annot


def MakeNoNa_RLavg3(fname,thresh = 0.5):
    file = pd.read_excel(fname, header = 0, index_col = 0)
    inum = file.index.to_frame().reset_index()
    jnum = file.columns.to_frame().reset_index()
    extend = fname.split('/')[1]
    param = extend.split('_')  #Macrophage_Endothelial_RLavg3mean_scaled
    lig = param[0]
    rec = param[1]
    for col in range(file.shape[1]):
        for row in range(file.shape[0]):
            if file.iloc[row,col] <= thresh:
                file.iloc[row,col] = 0 #p.nan
    file.to_excel('intermediate_choroid/'+lig+'_'+rec+'_'+'RLavg3mean_scaled_nona.xlsx')
    
    return True



def Cluster(fname,plot, cbar_ax = None,ax = None,fig = None,thresh = 0.5, minDist = 50, minSizeClust = 50, minSizeBounds = 50, getTop = None,drawBounds = True, findMax = True):
    extend = fname.split('/')[1]
    param = extend.split('_')  #Macrophage_Endothelial_RLavg3mean_scaled
    lig = param[0]
    rec = param[1]
    path = 'out_choroid/heatmaps/cluster'+str(thresh)+'/'
    makePath(path)
    if os.path.isfile(fname):
        file = pd.read_excel(fname, header = 0, index_col = 0)
    else:
        MakeNoNa_RLavg3('intermediate_choroid/'+lig+'_'+rec+'_RLavg3mean_scaled.xlsx',thresh = thresh)
        file = pd.read_excel(fname, header = 0, index_col = 0)
        
    
    try:
        b = pd.read_hdf('2D_Dict.h5', key = lig+'_'+rec+'_minDist_'+str(minDist)+'_minSizeClust_'+str(minSizeClust))
    except:
        b = MakeGroups2D(fname,minDist = minDist, minSizeClust = minSizeClust)
        pass
    c = pd.read_excel('intermediate_choroid/maxvals_topborders_DBSCAN.xlsx', index_col = 0)
    c.to_excel('intermediate_choroid/'+lig+'_'+rec+'_maxvals_topborders_DBSCAN.xlsx')
    if plot == True:
        file = pd.read_excel('intermediate_choroid/'+lig+'_'+rec+'_RLavg3mean_scaled_nona.xlsx',header = 0, index_col = 0)
        file = file.replace(0,np.nan)
        colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
        nrows = c.shape[0]

        possible_colors = list(colors.keys())
        
        #sns.set(rc={'figure.figsize':(26.7,23.27)}) #figure Size
        sns.set(font_scale = 0.05)
        ligand_order = file.index.unique()
        receptor_order = file.columns.unique()
        gEC_choroidLigands = ['ADIPOQ', 'ANGPT1', 'ANGPT2', 'ANGPTL3', 'BMP2', 'BMP7', 'C3', 'C4B', 'DLL1',
                              'FN1', 'IL13', 'IL4', 'JAG1', 'JAG2', 'LIF', 'OSM', 'S1P', 'SEMA3A', 'SEMA3C',
                              'SEMA3E', 'SEMA4A', 'SLIT2', 'TGFB1', 'TGFB2', 'TGFB3', 'VEGFA'] # 26

        ligands_1777 = np.loadtxt('geneLists/ligands_1777.txt',dtype =str).tolist() #"genes of interest"
        igands_44 = np.loadtxt('geneLists/ligands_44.txt', dtype=str).tolist() #"known regulators"

        interesting_ligands =  gEC_choroidLigands # "known regulators" from geneslist in DECNEO # ligands_44 if panglao, gEC_choroidLigands if choroid
        interesting_ligands_in = []
        interesting_ligands_index = []
        for x in interesting_ligands:
            try:
                for z in range(0,len(ligand_order)):
                    if x == ligand_order[z]:
                        interesting_ligands_index.append(z)
                        interesting_ligands_in.append(x)
            except:
                pass
    
        plt.figure(figsize = (15,16))
        if ax == None:
            ax1 = sns.heatmap(file,
                         cmap = sns.color_palette("YlGnBu", as_cmap=True),
                         linewidths = 0.05,linecolor = None,
                         xticklabels = True, yticklabels = True,
                         square = True,cbar = False, cbar_ax = cbar_ax,vmin = 0,vmax = 1)
        else:
            ax1 = sns.heatmap(file,
                    cmap = sns.color_palette("YlGnBu", as_cmap=True),
                    linewidths = 0.05,linecolor = None,
                    xticklabels = False, yticklabels = False,
                    square = False,cbar = False,cbar_ax = cbar_ax, vmin = 0,vmax = 1, ax = ax)
            
        fname_height = np.shape(file)[0]
        fname_width = np.shape(file)[1]
        #aspect equation: height = num * width

        #ax1.set_aspect(aspect = fname_height/fname_width)
        #ax1.xaxis.tick_top()  # x axis on top
        #ax1.xaxis.set_label_position('top')

        #plt.setp(ax1.get_xticklabels(), rotation=30, horizontalalignment='right')
        
        if drawBounds == True:
            try:
                c['rec_min'] = c['min_recs']
                c['rec_max'] = c['max_recs']
                c['lig_min'] = c['min_ligs']
                c['lig_max'] = c['max_ligs']
            except:
                pass
            for nrow in range(nrows):
                rec_min = c.iloc[nrow,].rec_min
                rec_max = c.iloc[nrow,].rec_max
                lig_min = c.iloc[nrow,].lig_min
                lig_max = c.iloc[nrow,].lig_max
                label_x,label_y = (rec_min+rec_max)/2 , (lig_min+lig_max)/2
                
                
                ax1.plot([rec_min,rec_min], [lig_min, lig_max], marker = 'o',ls = '--',c = possible_colors[nrow])
                ax1.plot([rec_max,rec_max], [lig_min,lig_max], marker = 'o',ls = '--',c = possible_colors[nrow])
                ax1.plot([rec_min,rec_max],[lig_min, lig_min], marker = 'o',ls = '--',c = possible_colors[nrow])
                ax1.plot([rec_min, rec_max], [lig_max,lig_max], marker = 'o',ls = '--',c = possible_colors[nrow])
                #ax1.text(label_x,label_y,c.index[nrow],c = possible_colors[nrow], size = 30,path_effects = [path_effects.Stroke(linewidth = 3, foreground = 'white'), path_effects.Normal()])

        #for x in range(0,len(interesting_ligands_in)-1):
        #    an = ax.annotate(interesting_ligands_in[x],xy = (10,interesting_ligands_index[x] + 0.5),xycoords = 'data',textcoords = 'data',
             #                xytext=(-150,interesting_ligands_index[x]+4), fontsize = 30,
              #               arrowprops=dict(arrowstyle="fancy", color = 'black', alpha = 0.2))
        #plt.tight_layout()
        #ax1.figure.savefig(path+'cluster_'+lig+'_'+rec+'.png',dpi = 1000)
    return ax1

