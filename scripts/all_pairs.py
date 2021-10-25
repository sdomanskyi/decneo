import os
import sys
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import scipy.cluster.hierarchy
from matplotlib import cm
from decneo.commonFunctions import read, write
import multiprocessing

cwd = '/mnt/gs18/scratch/users/paterno1/otherCellTypes_choroid/LRpairs/'

celltypes = ['Endothelial', 'Pericyte', 'Fibroblast', 'Macrophage', 'SMC']
genetypes = ['ligands', 'receptors']

def loadRamiowskiLRPairs(dir):
    df = pd.read_excel(dir + 'SupplementaryData2Ramilowski.xlsx', sheet_name='All.Pairs', index_col=False)
    df = df[['Ligand.ApprovedSymbol', 'Receptor.ApprovedSymbol', 'Pair.Evidence']]
    df = df.loc[df['Pair.Evidence'] != 'EXCLUDED']
    df = df.loc[df['Pair.Evidence'] != 'EXCLUDED not receptor']
    df = df.loc[df['Pair.Evidence'] != 'EXCLUDED not ligand']
    df = df.drop(columns=['Pair.Evidence'])
    df.columns = ['ligand', 'receptor']
    df = df.reset_index(drop=True)
    return pd.MultiIndex.from_arrays([df['ligand'].values, df['receptor'].values], names=['ligand', 'receptor'])

def doLR(geneL, geneR, cutoff, hcutoff, dataDict, LR = None, suffix = 'data', makePlot = True, saveData = True):
    
    name = '%s-%s' % (geneL, geneR)
    
    grays = plt.cm.gray_r
    
    resDict = dict()
    for genetype in genetypes:
        resDict.update({genetype: dict()})    
        for celltype in celltypes:   
            df_peaks = dataDict[hcutoff][genetype][celltype]
            gene = geneL if genetype=='ligands' else geneR
            se_peak = df_peaks[gene] if gene in df_peaks.columns else pd.Series(dtype=float)
            resDict[genetype].update({celltype: se_peak[se_peak>=cutoff]})
            
            
    if makePlot:
        fig, ax = plt.subplots(figsize=[7,7])
    cx = np.array([0.,0.25,0.5,0.75,1.])*1.35
    gy = [0.2, 0.8]
    colorg = ['blue', 'green']
    sh = [-1, 1]
    de = 0.035
    celltypesO = ['SMC', 'Pericyte', 'Endothelial', 'Fibroblast', 'Macrophage']

    if makePlot:
        ax.text(0, 0.9, name, ha='center', va='center', fontsize=20)
    
    gg = []
    for ig1, genetype1 in enumerate(genetypes):
        if makePlot:
            ax.text(0.5*1.35, gy[ig1]+0.05*sh[ig1], genetype1, ha='center', va='center', fontsize=20)

        for ic1, celltype1 in enumerate(celltypesO):
            if makePlot:
                ax.text(cx[ic1], gy[ig1]+0.00*sh[ig1], celltype1, ha='center', va='center', fontsize=15)

            t1 = resDict[genetype1][celltype1]
            group1 = 0
            h1 = 1.
            g1 = t1
            temp1 = cx[ic1] + (-1/2 + group1 + 0.5)*de, gy[ig1]-0.05*sh[ig1]

            if makePlot:
                mec = 'k'
                ax.plot(*temp1, 'o', ms=0.9*len(g1)/2, color=colorg[ig1], mec=mec, mew=1.0)

            ggc = []
            ig2, genetype2 = 1, 'receptors'
            if genetype2!=genetype1:
                for ic2, celltype2 in enumerate(celltypesO):
                    t2 = resDict[genetype2][celltype2]
                    group2 = 0
                    temp2 = cx[ic2] + (1/2-group2-0.5)*de, gy[ig2]-0.05*sh[ig2]

                    c = pd.MultiIndex.from_tuples([(a, b) for a in g1.index for b in t2.index], names=['ligand', 'receptor'])
                    comP = LR.intersection(c)
                    com = comP.shape[0]

                    if makePlot:
                        if com>0:
                            alpha = 1. if com >= 15 else max(com/15., 0.2)
                            ax.annotate("", xy=(temp2[0], temp2[1]), xycoords='data', xytext=(temp1[0], temp1[1]), textcoords='data', 
                                        arrowprops=dict(facecolor=grays(alpha), edgecolor=grays(alpha), shrink=0.04, headwidth=7, width=com/5, alpha=1.), zorder=-1)

                    ggc.append(', '.join(list(comP.get_level_values(0) + '-' + comP.get_level_values(1))) if com>0 else None)

            templ = [genetype1, celltype1, ', '.join(list(g1.index.values))]
            templ.extend(ggc)
            gg.append(templ)

    if makePlot:
        plt.axis('off')
        fig.tight_layout()

        plt.savefig('%s %s.png' % (name, suffix), dpi=300)
        exportEmf(cwd, '%s %s' % (name, suffix))
    
    dfx = pd.DataFrame(gg).drop(0, axis=1)
    dfx.columns = ['cell type', 'ligands', 'SMC pairs', 'Pericyte pairs', 'Endothelial pairs', 'Fibroblast pairs', 'Macrophage pairs']
    recc = dfx.iloc[5:]['ligands'].values
    dfx = dfx.iloc[:5]
    dfx.insert(2, 'receptors', recc)

    return dfx

def doLRw(params):
    args, kwargs = params
    return doLR(*args, **kwargs).set_index('cell type', drop=True).stack().rename((args[0], args[1]))

if __name__ == '__main__':

    LR = loadRamiowskiLRPairs(cwd)
    dataDictVoigt = read(cwd + 'dataDictVoigtv2_10_14_2021')  
    dataDictPanglao = read(cwd + 'dataDictPanglao')

    print('LR pairs from Ramilowski et al. 2015:\n', LR, flush=True)

    if True:
        for hcutoff in ['0.0']:
            allDataDicts = {'choroid': dataDictVoigt}
            #allDataDicts = {'choroid': dataDictVoigt, 'mouse': dataDictPanglao['Mus musculus'], 'human': dataDictPanglao['Homo sapiens']}
            for key in allDataDicts.keys():
                tligands = pd.concat([allDataDicts[key][hcutoff]['ligands'][celltype].index.to_series() for celltype in celltypes]).drop_duplicates().sort_values().index
                treceptors = pd.concat([allDataDicts[key][hcutoff]['receptors'][celltype].index.to_series() for celltype in celltypes]).drop_duplicates().sort_values().index
                tpairs = pd.Series(index=pd.MultiIndex.from_product([tligands, treceptors], names=['ligand', 'receptor']), data=0)
                tpairs = tpairs.index

                nCPUs = 100
                print('hcutoff:', hcutoff, ', data:', key, 'nCPUs:', nCPUs, ', pairs:', len(tpairs))

                pool = multiprocessing.Pool(processes=nCPUs)
                df = pd.DataFrame(pool.map(doLRw, [((pair[0], pair[1], 0.3, hcutoff, allDataDicts[key]), dict(LR=LR, makePlot=False, saveData=False)) for pair in tpairs.values]))
                pool.close()
                pool.join()

                df.to_hdf('v3_allpairs_%s_%s.h5' % (hcutoff, key), key='df', mode='a', complevel=4, complib='zlib')
                del df

    if True:
        for hcutoff in ['0.0']:
            allDataDicts = {'choroid': dataDictVoigt}
            #allDataDicts = {'choroid': dataDictVoigt, 'mouse': dataDictPanglao['Mus musculus'], 'human': dataDictPanglao['Homo sapiens']}
            for key in allDataDicts.keys():
                print('hcutoff:', hcutoff, ', data:', key)
                df_all = pd.read_hdf(cwd + 'v3_allpairs_%s_%s.h5' % (hcutoff, key), key='df')
                for celltype1 in celltypes:
                    for celltype2 in celltypes:
                        # celltype1 has ligands, celltype2 has receptors
                        df = df_all.xs([(celltype1, celltype2 + ' pairs')], axis=1).fillna('')
                        tempName = cwd + 'v3_allpairs_%s_%s_%s_%s.h5' % (hcutoff, key, celltype1, celltype2)
                        df.to_hdf(tempName + '.h5', key='df', mode='a', complevel=4, complib='zlib')   
                        df = pd.read_hdf(tempName + '.h5', key='df')
                        l = len(df)
                        df.index = pd.MultiIndex.from_tuples(df.index.values, names=['ligand', 'receptor'])
                        df.columns = ['pairs']
                        df['count'] = df['pairs'].apply(lambda v: len(set(v.split(', ')) - {''}))
                        df = df.sort_values(by='count', ascending=False)
                        df['ram'] = ~(pd.Series(index=LR, data=0).reindex(df.index).isna())
                        df = df.loc[df['count']>0]
                        df['count-1'] = df['count'] - df['ram'].astype(int)
                        df = df[['ram', 'count', 'count-1', 'pairs']]
                        df.reset_index().to_excel(tempName + '.xlsx', sheet_name=str(l) + ' pairs tested', index=False)

                        df = df['count-1']
                        for mode in ['ligand', 'receptor']:
                            df.groupby(level=mode).max().sort_values(ascending=False).to_excel(tempName[:-3] + '_%s_max.xlsx' % mode)
