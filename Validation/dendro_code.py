import os
import numpy as np
import pandas as pd
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from adjustText import adjust_text
from sklearn.decomposition import PCA
plt.ioff()


def dendro(df_corr, 
           #selGenes, 
           
           #celltype, 
           genesSubset,
           fig_file = None,
           inhLoc = 0, 
           metric = 'euclidean', 
           linkageMethod = 'ward',
           pca = None,
           n_clusters=5):

    fig, axs = plt.subplots(1, 1, figsize=(15, 8))

    #for i, species in enumerate(['Homo sapiens', 'Mus musculus']):

    #print('Reading corelations of %s data' % species, flush=True)
    #df_corr = pd.read_hdf(corrFile, key=species)
    #df_corr = df_corr[df_corr.columns.intersection(np.unique(selGenes))]
    print(df_corr.shape, flush=True)

    ax = axs#[0]

    origLineWidth = matplotlib.rcParams['lines.linewidth']
    matplotlib.rcParams['lines.linewidth'] = 0.75
    
    if pca == None:
        Z = hierarchy.linkage(df_corr.values.T, metric=metric, method=linkageMethod, optimal_ordering=True)
    else:
         Z = hierarchy.linkage(PCA(n_components= pca).fit_transform(df_corr.values.T), metric=metric, method=linkageMethod, optimal_ordering=True)
    
    cmap = cm.gist_ncar(np.linspace(0, 0.5, n_clusters + 1)) #gist_ncar #nipy_spectral #hsv
    hierarchy.set_link_color_palette([matplotlib.colors.rgb2hex(rgb[:3]) for rgb in cmap])
    D = hierarchy.dendrogram(Z, ax=ax, color_threshold = (Z[-n_clusters,2] + Z[-n_clusters+1,2]) / 2, above_threshold_color='lightgrey', orientation='top')
    

    reindexed = pd.Index(df_corr.columns[D['leaves']]).reindex(pd.Index(genesSubset).intersection(df_corr.columns))
    
    hierarchy.set_link_color_palette(None)

    matplotlib.rcParams['lines.linewidth'] = origLineWidth
    
    genes = reindexed[0][reindexed[1] > -1].values
    locations = reindexed[1][reindexed[1] > -1]

    tickLabelsColors = np.array(['black']*len(df_corr.columns), dtype=np.dtype('U20'))
    xtickslabels = np.array(['']*len(df_corr.columns), dtype=np.dtype('U20'))
    for gene, location in zip(genes, locations):
        xtickslabels[location] = gene
        tickLabelsColors[location] = 'green' if (gene in genesSubset[:inhLoc]) else 'red'

    ax.set_xticklabels(xtickslabels, fontsize=3.5)

    for xtick, color in zip(ax.get_xticklabels(), tickLabelsColors):
        xtick.set_color(color)

    texts = []
    origPos = []
    for xpos, xtext, color in zip(ax.get_xticks(), xtickslabels, tickLabelsColors):
        if xtext != '':
            texts.append(ax.text(xpos, -2., xtext, fontsize=6, rotation=90, va='top', ha='center', color=color))
            origPos.append(xpos)

    adjust_text(texts, va='top', ha='center', autoalign='x', lim=200, only_move={'text':'x'})

    v = 0.04 * ax.get_ylim()[1]

    for text, opos in zip(texts, origPos):
        text._y = -v
        ax.plot([text._x, opos], [text._y, 0.], color=text._color, lw=0.5, clip_on=False)

    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])
    #ax.set_ylabel(species)

    #fig.suptitle('Celltype: %s' % celltype, fontsize=12)
    #if fig != None:
    #    fig.savefig(fig_file, dpi=600)
    fig.savefig("tmp")
    
    reindexed = pd.Index(df_corr.columns[D['leaves']]).reindex(df_corr.columns)
    reindexed = pd.DataFrame(reindexed[1],index = reindexed[0], columns = ["Dendrogram"])

    return reindexed