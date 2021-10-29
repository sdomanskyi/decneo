import numpy as np
import copy
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import sklearn.cluster
from sklearn.cluster import DBSCAN
from decneo.commonFunctions import read, write,adjustTexts1D
import os
from adjustText import adjust_text

def makePath(path):
    if not os.path.exists(path):
        os.makedirs(path)
        
    return True
def countUnique(s):
    return len([0 for sub in s.values[0].split(', ') if (not s.name[0] in sub.split('-')) and (not s.name[1] in sub.split('-'))])

def findGroupsDBSCAN(dfM, eps, min_samples, figwidth = None, ax = None, plot = True, findMaxVal = False, fpath = None, getTop = None):
    dfM = dfM #[::-1]
    inum = dfM.index.to_frame().reset_index()
    jnum = dfM.columns.to_frame().reset_index()
    data = dfM.replace(0, np.nan).reset_index(drop=True).T.reset_index(drop=True).T.stack().index.to_frame().values
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(data)
    uclusters = pd.Series(db.labels_).value_counts()
    dfL = pd.Series(index=pd.MultiIndex.from_arrays(data.T), data=db.labels_).unstack(-1).reindex(inum.index, axis=0).reindex(jnum.index, axis=1)
    
    if plot == True:
        if ax is None:
            fig, ax = plt.subplots(figsize=(figwidth, figwidth))
        
        ligand_order = dfM.index.unique()
        receptor_order = dfM.columns.unique()
        gEC_choroidLigands = ['ADIPOQ', 'ANGPT1', 'ANGPT2', 'ANGPTL3', 'BMP2', 'BMP7', 'C3', 'C4B', 'DLL1',
                              'FN1', 'IL13', 'IL4', 'JAG1', 'JAG2', 'LIF', 'OSM', 'S1P', 'SEMA3A', 'SEMA3C',
                              'SEMA3E', 'SEMA4A', 'SLIT2', 'TGFB1', 'TGFB2', 'TGFB3', 'VEGFA'] # 26


        interesting_ligands =  gEC_choroidLigands # "known regulators" from geneslist in DECNEO # ligands_44 if panglao, gEC_choroidLigands if choroid
        receptorsListHugo_2555 = np.loadtxt('geneLists/receptorsListHugo_2555.txt', dtype=str).tolist()
        receptorsListHugo = np.loadtxt('geneLists/receptors_Human_HUGO.txt', dtype=str).tolist()
        receptorsListHGNC = np.loadtxt('geneLists/HGNC_receptors.txt', dtype=str).tolist()
        gEC22 = ['KDR','FLT1','FLT4','NRP1','NRP2','FGFR1','FGFR2','FGFR3','CXCR2','ROBO1',
         'ROBO4','ENG','PDGFRA','PDGFRB','TEK','KIT','MET','CLEC14A', # stimulators
         'CD36','CD47','VLDLR','PLXND1']
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
            
        interesting_receptors = gEC22
        interesting_receptors_in = []
        interesting_receptors_index = []
        
        for x in interesting_receptors:
            try:
                for z in range(0,len(receptor_order)):
                    if x == receptor_order[z]:
                        interesting_receptors_index.append(z)
                        interesting_receptors_in.append(x)
            except:
                pass
        
        cmap = copy.copy(plt.cm.jet); cmap.set_bad('white'); cmap.set_under('grey')
        pdata = dfL.values #[::-1]
        ax.pcolormesh(pdata, cmap=cmap, vmin=0)
        ax.tick_params(axis="x", bottom=False, labelbottom = False, top=True, labeltop=True)
        #ax.set_xticklabels([]); ax.set_xticks([]); ax.set_yticklabels([]); ax.set_yticks([]);
        # print((~dfL.isna()).sum().sum()/(dfL.shape[0]*dfL.shape[1]))

        resDict = dict()
        dftemp = dfM.reset_index(drop=True).T.reset_index(drop=True).T
        min_ligs, max_ligs, min_recs, max_recs, max_vals = dict(),dict(),dict(),dict(), dict()
        
        for ucluster in uclusters.index:
            if ucluster>-1:
                wh = np.where((dfL.fillna(-2).astype(int)==ucluster))
                min_ligs[ucluster], max_ligs[ucluster] = wh[0].min(), wh[0].max() + 1
                min_recs[ucluster], max_recs[ucluster] = wh[1].min(), wh[1].max() + 1
                max_vals[ucluster] = np.max(np.max(dfM.iloc[min_ligs[ucluster]:max_ligs[ucluster],min_recs[ucluster]:max_recs[ucluster]]))
                cmean = dftemp.values[wh].mean()
                boxmean = dftemp.iloc[wh[0].min(): wh[0].max()+1, wh[1].min(): wh[1].max()+1].mean().mean()
                boxsize = (wh[0].max() - wh[0].min()) * (wh[1].max() - wh[1].min())
                resDict.update({ucluster: {'intensity': cmean, 'boxintensity': boxmean, 'size': len(wh[0]), 'boxsize': boxsize, 
                                       'ligands': str(np.sort(inum.loc[wh[0]]['index'].unique()).tolist()).strip('][').replace("'",''), 
                                       'boxligands': str(np.sort(inum.iloc[wh[0].min(): wh[0].max()+1]['index'].unique()).tolist()).strip('][').replace("'",''), 
                                       'receptors': str(np.sort(jnum.loc[wh[1]]['index'].unique()).tolist()).strip('][').replace("'",''),
                                       'boxreceptors': str(np.sort(jnum.iloc[wh[1].min(): wh[1].max()+1]['index'].unique()).tolist()).strip('][').replace("'",'')}})
                
        out = pd.DataFrame(resDict).T.sort_values(by='boxintensity', ascending=False)
        #print(out)
        df = pd.concat([pd.Series(min_recs),pd.Series(max_recs),pd.Series(min_ligs),pd.Series(max_ligs),pd.Series(max_vals)],axis = 1)
        df.columns = ['min_recs','max_recs','min_ligs','max_ligs','max_vals']
        df = df.sort_values(by = 'max_vals', ascending = False)
        if findMaxVal == True:
            row_max_val = df.sort_values(by = 'max_vals', ascending = False).index[0]
            print(row_max_val)
            rec_min,rec_max,lig_min,lig_max = min_recs[row_max_val],max_recs[row_max_val],min_ligs[row_max_val],max_ligs[row_max_val]
            ax.plot([rec_min,rec_min], [lig_min, lig_max], marker = 'o',ls = '--', c = 'black')
            ax.plot([rec_max,rec_max], [lig_min,lig_max], marker = 'o',ls = '--', c = 'black')
            ax.plot([rec_min,rec_max],[lig_min, lig_min], marker = 'o',ls = '--', c = 'black')
            ax.plot([rec_min, rec_max], [lig_max,lig_max], marker = 'o',ls = '--', c = 'black')
            #out = pd.DataFrame(out.loc[row_max_val,:]).T
        elif findMaxVal == False and getTop == None:
            for nrow in range(np.shape(df)[0]):
                rec_min,rec_max,lig_min,lig_max = min_recs[nrow],max_recs[nrow],min_ligs[nrow],max_ligs[nrow]
                ax.plot([rec_min,rec_min], [lig_min, lig_max], marker = 'o',ls = '--', c = 'black')
                ax.plot([rec_max,rec_max], [lig_min,lig_max], marker = 'o',ls = '--', c = 'black')
                ax.plot([rec_min,rec_max],[lig_min, lig_min], marker = 'o',ls = '--', c = 'black')
                ax.plot([rec_min, rec_max], [lig_max,lig_max], marker = 'o',ls = '--', c = 'black')
        
        if getTop != None:
            row_max_vals = df.sort_values(by = 'max_vals', ascending = False).index[0:getTop].tolist()
            df_maxvals_top = df.sort_values(by = 'max_vals', ascending = False).iloc[0:getTop]
            df_maxvals_top.to_excel('intermediate_choroid/maxvals_topborders_DBSCAN.xlsx')
            for nrow in row_max_vals:
                rec_min,rec_max,lig_min,lig_max = min_recs[nrow],max_recs[nrow],min_ligs[nrow],max_ligs[nrow]
                ax.plot([rec_min,rec_min], [lig_min, lig_max], marker = 'o',ls = '--', c = 'black')
                ax.plot([rec_max,rec_max], [lig_min,lig_max], marker = 'o',ls = '--', c = 'black')
                ax.plot([rec_min,rec_max],[lig_min, lig_min], marker = 'o',ls = '--', c = 'black')
                ax.plot([rec_min, rec_max], [lig_max,lig_max], marker = 'o',ls = '--', c = 'black')
                ax.text(y = (lig_min + lig_max)/2, x = (rec_min+rec_max)/2, s = nrow, size = 30, path_effects = [path_effects.Stroke(linewidth = 3, foreground = 'white'), path_effects.Normal()])
            print(df)
            top = df.iloc[0:getTop,].index
            out = out.loc[top,:]
            print(out)
        if fpath != None:
            temp_str = fpath.split('/')
            folder = '/'.join(temp_str[0:len(temp_str)-1])+'/'
            makePath(folder)
    else:
        out = pd.DataFrame(resDict).T.sort_values(by='boxintensity', ascending=False)
            
    return out

def Group_DBSCAN(fname,eps, min_samples,fig = None, ax = None, figwidth = 200, plot = True, findMaxVal = True,getTop = None):
    lig = fname.split('/')[1].split('_')[0]
    rec = fname.split('/')[1].split('_')[1]
    fpath = 'out_choroid/heatmaps/DBSCAN/'+lig+'_'+rec+'_RLavg3.png' #'intermediate_choroid/'+lig_name+'_'+rec_name+'_RLavg3mean_scaled_nona.xlsx'

    temp_str = fpath.split('/')
    folder = '/'.join(temp_str[0:len(temp_str)-1])+'/'
    makePath(folder)
    dfM = pd.read_excel(fname,index_col = 0, header = 0)
    M = dfM.values
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(15, 16))
    else:
        ax = ax
        fig = fig
    cmap = copy.copy(plt.cm.jet_r); cmap.set_bad('white')
    cM = M.copy()
    #params = dict(cmap=cmap, vmin=0, vmax=50, aspect = 'equal',interpolation='None', extent=(-0.5, cM.shape[0] - 0.5, cM.shape[1] - 0.5, -0.5)) #aspect = 'auto'
    params = dict(cmap=cmap, vmin=0, vmax=50, aspect = 'auto',interpolation='None', extent=(0, cM.shape[0], cM.shape[1],0))
    cax = ax
    cax.imshow(np.ma.array(cM, mask=(cM==0)), **params)
    
    
    rec_labels = [rec for rec in dfM.columns]
    lig_labels = [lig for lig in dfM.index]
    lig_labels = lig_labels[::-1]
    
     
    
    cax.set_xticks(np.arange(0, len(rec_labels), 1.0)); cax.set_xlim(0, len(rec_labels)); cax.set_xticklabels(rec_labels, fontsize = 1, rotation=90);  cax.set_yticks(np.arange(0, len(lig_labels),1.0)); cax.set_ylim(0, len(lig_labels)); cax.set_yticklabels(lig_labels, fontsize = 1);    
    
    out = findGroupsDBSCAN(dfM, eps = eps, min_samples = min_samples, figwidth = figwidth, ax = cax, plot = plot, findMaxVal = findMaxVal, fpath = fpath, getTop = getTop);
    out.to_excel(folder+lig+'_'+rec+'_RLavg3.xlsx')
    return True

