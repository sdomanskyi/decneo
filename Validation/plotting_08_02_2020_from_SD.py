import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
import matplotlib
from scipy.cluster import hierarchy
from matplotlib import cm
from adjustText import adjust_text
import scipy
import matplotlib.patheffects as path_effects
from scipy.spatial.distance import squareform
from sklearn.metrics.pairwise import euclidean_distances
import os
def plotting(df, 
            stimulators,
            inhibitors,
            bar_df,
            saveDir,
            metric = 'euclidean', 
            linkageMethod = 'ward', 
            n_clusters = 10, 
            adjustText = True,
            majorMetric = "Correlation",
            suffix = ""):

    '''df has all genes in rows, and receptors in columns, values are gene expression correlation'''
                                                
    def addDendro(fig, dataGenes, M, coords, linewidth=0.25, adjustText = adjustText):

        genesSubset = list(stimulators) + list(inhibitors)

        ax = fig.add_axes(coords, frame_on=False)
        Z = hierarchy.linkage(np.nan_to_num(M, nan=max(M)), method=linkageMethod, optimal_ordering=True)

        origLineWidth = matplotlib.rcParams['lines.linewidth']
        matplotlib.rcParams['lines.linewidth'] = linewidth
        cmap = cm.gist_ncar(np.linspace(0, 0.5, n_clusters + 1))
        hierarchy.set_link_color_palette([matplotlib.colors.rgb2hex(rgb[:3]) for rgb in cmap])

        D = hierarchy.dendrogram(Z, ax=ax, color_threshold = (Z[-n_clusters,2] + Z[-n_clusters+1,2]) / 2, above_threshold_color='k', orientation='top')
        hierarchy.set_link_color_palette(None)

        matplotlib.rcParams['lines.linewidth'] = origLineWidth

        reindexed = pd.Index(dataGenes[D['leaves']]).reindex(pd.Index(genesSubset).intersection(dataGenes))
        genes = reindexed[0][reindexed[1] > -1].values
        locations = reindexed[1][reindexed[1] > -1]

        if True:
            tickLabelsColors = np.array(['navy']*len(dataGenes), dtype=np.dtype('U20'))
            xtickslabels = np.array(['']*len(dataGenes), dtype=np.dtype('U20'))
            for gene, location in zip(genes, locations):
                xtickslabels[location] = gene
                tickLabelsColors[location] = 'green' if (gene in stimulators) else 'red'

            ax.set_xticklabels(xtickslabels, fontsize=4)
            ax.tick_params(axis='y', labelsize=4, width=0.25, length=1)
                    
            ax.set_yticklabels([])
            ax.set_yticks([])

            for xtick, color in zip(ax.get_xticklabels(), tickLabelsColors):
                xtick.set_color(color)

            texts = []
            origPos = []
            for xpos, xtext, color in zip(ax.get_xticks(), xtickslabels, tickLabelsColors):
                if xtext != '':
                    texts.append(ax.text(xpos, -2., xtext, fontsize=6, rotation=90, va='top', ha='center', color=color))
                    origPos.append(xpos)

            ticks_x = []
            ticks_y = []
            vdistance = -0.01 * ax.get_ylim()[1]
            for tick in ax.get_xticks():
                ticks_x.extend([tick, tick, None])
                ticks_y.extend([0, vdistance, None])
            ax.plot(ticks_x, ticks_y, color='k', lw=0.4, clip_on=False)
            ax.set_xticklabels([])
                        

            if adjustText:
                adjust_text(texts, va='top', ha='center', autoalign='x', lim=400, only_move={'text':'x'})

            v = 0.04 * ax.get_ylim()[1]
            for text, opos in zip(texts, origPos):
                text._y = -v
                ax.plot([text._x, opos], [text._y, 0.], color=text._color, lw=0.5, clip_on=False)

        if True:
            clusters = scipy.cluster.hierarchy.fcluster(Z, t=n_clusters, criterion='maxclust')[D['leaves']] - 1
            clusterBoundaries = (np.where(clusters - np.roll(clusters, 1) != 0)[0]/ len(D['leaves'])) * ax.get_xlim()[1]
            clusterBoundaries = np.append(clusterBoundaries, ax.get_xlim()[1])

            clusterCenters = clusterBoundaries[:-1] + ((clusterBoundaries - np.roll(clusterBoundaries, 1))/2.)[1:]
            vposition = (Z[-n_clusters,2] + Z[-n_clusters+1,2]) / 5

            for cluster, position in zip(np.unique(clusters), clusterCenters):
                ltext = ax.text(position, vposition, '#%s' % cluster, fontsize=7, color='white', va='center', ha='center')
                ltext.set_path_effects([path_effects.Stroke(linewidth=1., foreground='k'), path_effects.Normal()])

        return {'order': D['leaves'], 
                'M': squareform(M)[:, D['leaves']][D['leaves'], :], 
                'genes': genes, 
                'allGenes': dataGenes[D['leaves']], 
                'locations': locations,
                'tickLabelsColors': tickLabelsColors,
                'xtickslabels': xtickslabels,
                'clusters': clusters,
                'clusterBoundaries': clusterBoundaries / 10.,
                'clusterCenters': clusterCenters / 10.}

    def addHeatmap(fig, dataArgs, coords, adjustText = adjustText):

        M = dataArgs['M'] 
        order = dataArgs['order'] 
        genes = dataArgs['genes'] 
        locations = dataArgs['locations'] 
        tickLabelsColors =  dataArgs['tickLabelsColors'] 
        tickslabels =  dataArgs['xtickslabels'] 
        clusters =  dataArgs['clusters'] 
        clusterBoundaries =  dataArgs['clusterBoundaries']
        clusterCenters =  dataArgs['clusterCenters']

        ax = fig.add_axes(coords, frame_on=False)

        masked_M = np.ma.array(M, mask=np.isnan(M))
        cmap = plt.cm.Greens_r
        cmap.set_bad('red')

        im = ax.imshow(masked_M, cmap=cmap, aspect='auto', interpolation='None', extent=(-0.5, M.shape[0] - 0.5, M.shape[1] - 0.5, -0.5))

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        # Selected x tick labels
        if True:
            ax.set_xticks(range(len(tickslabels)))
            ax.set_xticklabels(tickslabels, fontsize=4)
            for xtick, color in zip(ax.get_xticklabels(), tickLabelsColors):
                xtick.set_color(color)

            texts = []
            origPos = []
            for xpos, xtext, color in zip(ax.get_xticks(), tickslabels, tickLabelsColors):
                if xtext != '':
                    texts.append(ax.text(xpos, 1.01*ax.get_ylim()[0], xtext, fontsize=6, rotation=90, va='top', ha='center', color=color))
                    origPos.append(xpos)

            ax.set_xticklabels([])
            ax.set_xticks([])

            if adjustText:
                adjust_text(texts, va='top', ha='center', autoalign='x', lim=400, only_move={'text':'x'})

            v = ax.get_ylim()[0]
            for text, opos in zip(texts, origPos):
                text._y = 1.01 * v
                ax.plot([text._x, opos], [text._y, v], color=text._color, lw=0.5, clip_on=False)
        
        # Selected y tick labels
        if True:
            ax.set_yticks(range(len(tickslabels)))
            ax.set_yticklabels(tickslabels, fontsize=4)
            for ytick, color in zip(ax.get_yticklabels(), tickLabelsColors):
                ytick.set_color(color)

            texts = []
            origPos = []
            for ypos, xtext, color in zip(ax.get_yticks(), tickslabels, tickLabelsColors):
                if xtext != '':
                    texts.append(ax.text(-0.01*ax.get_xlim()[1], ypos, xtext, fontsize=6, va='center', ha='right', color=color))
                    origPos.append(ypos)

            ax.set_yticklabels([])
            ax.set_yticks([])

            if adjustText:
                adjust_text(texts, va='center', ha='right', autoalign='y', lim=400, only_move={'text':'y'})

            v = -0.01 * ax.get_xlim()[1]
            for text, opos in zip(texts, origPos):
                text._x = v
                ax.plot([0., text._x], [opos, text._y], color=text._color, lw=0.5, clip_on=False)

        # Clusters outline boxes
        if True:
            for cluster, position in zip(np.unique(clusters), clusterCenters):
                ltext = ax.text(position, position, '#%s' % cluster, fontsize=7, color='white', va='center', ha='center')
                ltext.set_path_effects([path_effects.Stroke(linewidth=1., foreground='k'), path_effects.Normal()])

            clusterBoundaries -= 0.5

            for i in range(len(np.unique(clusters))):
                ax.plot([clusterBoundaries[i], clusterBoundaries[i+1], clusterBoundaries[i+1], clusterBoundaries[i], clusterBoundaries[i]], 
                        [clusterBoundaries[i], clusterBoundaries[i], clusterBoundaries[i+1], clusterBoundaries[i+1], clusterBoundaries[i]], 
                        '--', lw=0.75, color='k', clip_on=False)

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        # Colorbar
        if True:
            ax = fig.add_axes([0.85, 0.1, 0.025, 0.6], frame_on=False)

            ax.set_xticks([])
            ax.set_xticklabels([])
            ax.set_yticks([])
            ax.set_yticklabels([])

            clb = fig.colorbar(im, ax=ax, fraction=0.4, label='Eucl. dist. of gene expr. %s dist.' % majorMetric)
            clb.ax.tick_params(labelsize=6)

        return

    def addBar(fig, dataArgs, mode, coords):

        M = dataArgs['M'] 
        order = dataArgs['order'] 
        genes = dataArgs['genes'] 
        locations = dataArgs['locations'] 
        tickLabelsColors =  dataArgs['tickLabelsColors'] 
        tickslabels =  dataArgs['xtickslabels'] 
        clusters =  dataArgs['clusters'] 
        clusterBoundaries =  dataArgs['clusterBoundaries']
        clusterCenters =  dataArgs['clusterCenters']
        allGenes = dataArgs['allGenes'] 

        ax = fig.add_axes(coords, frame_on=True)
        ax.set_xlim([min(clusterBoundaries), max(clusterBoundaries)])

        if mode == 0:
            ylabel='Binomial\nP-Val'
            data =pd.read_hdf(bar_df, key='df')["BN"].reindex(allGenes).values

        elif mode == 1:
            ylabel='Upregulated\nP-Val'
            data =pd.read_hdf(bar_df, key='df')["P-Val"].reindex(allGenes).values

        elif mode == 2:
            ylabel='Fold\nChange'
            data =pd.read_hdf(bar_df, key='df')["LFC"].reindex(allGenes).values

        elif mode == 3:
            ylabel='Fraction'
            data =pd.read_hdf(bar_df, key='df')["Fraction"].reindex(allGenes).values

        elif mode == 4:
            clust_df = pd.DataFrame(clusters,index = allGenes,columns = ["Clust"])
            #clust_df = clust_df.sort_values("Order")
            #clust_df["Clust"] = cluster
            all_values = []
            for c in set(clust_df["Clust"]):
                cg = clust_df.loc[clust_df["Clust"] == c].index
                
                euclid_df = euclidean_distances(df[cg].fillna(0).T)
                euclid_df = pd.DataFrame(euclid_df,index= cg,columns=cg)
                euclid_values = euclid_df.mean().values
                all_values += list(euclid_values)
            data = [max(all_values) - x for x in all_values]
            ylabel = "Cluster\nCloseness"
            
        elif mode == 5:
            clust_df = pd.DataFrame(clusters,index = allGenes,columns = ["Clust"])
            #clust_df = clust_df.sort_values("Order")
            #clust_df["Clust"] = cluster
            all_values = []
            euclid_df = euclidean_distances(df.loc[allGenes,allGenes].fillna(0).T)
            euclid_df = pd.DataFrame(euclid_df,allGenes,allGenes)
            all_values = []
            for g in allGenes:
                neigh = dendro_dist.sort_values(g).head(21).tail(20).index
                if g == "KDR":
                    print(neigh)
                all_values.append(euclid_df.loc[g,neigh].mean())
                
               
            data = [max(all_values) - x for x in all_values]
            ylabel = "Neighborhood\nCloseness"
        
        elif mode == 6:
            ylabel='Angiogenesis\nLiterature'
            data =pd.read_hdf(bar_df, key='df')["angiogenesis"].reindex(allGenes).values


        elif mode == 7:
            ylabel='Endothelial\nLiterature'
            data =pd.read_hdf(bar_df, key='df')["endothelial"].reindex(allGenes).values

        elif mode == 8:
            ylabel='Conservation'
            data =pd.read_hdf(bar_df, key='df')["DD-Conservation"].reindex(allGenes).values
            
            
        elif mode == 9:
            ylabel='ED_Conservation'
            data =pd.read_hdf(bar_df, key='df')["ED-Conservation"].reindex(allGenes).values
         
        elif mode == 10:
            ylabel='CC_Conservation'
            data =pd.read_hdf(bar_df, key='df')["CC-Conservation"].reindex(allGenes).values
        
        
        elif mode == 11:
            ylabel='All 4 Window\nAvaerage'
            data =pd.read_hdf(bar_df, key='df')["All-4_WS21"].reindex(allGenes).values
            
        elif mode == 12:
            ylabel='Ind 3 Window\nAvaerage'
            data =pd.read_hdf(bar_df, key='df')["Independent-3_WS21"].reindex(allGenes).values


        ax.bar(range(len(clusters)), data, width=ax.get_xlim()[1]/len(clusters), color=tickLabelsColors)
        ax.set_xticks([])
        ax.set_xticklabels([])

        yticks = np.round(ax.get_ylim(), 1)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks)

        ax.tick_params(axis='y', labelsize=6, width=0.75, length=3)

        if True:
            ylim = ax.get_ylim()
            for i in range(1, len(np.unique(clusters))):
                ax.plot([clusterBoundaries[i] - 0.5]*2, [ylim[0], ylim[1]], '--', lw=0.5, color='k', clip_on=False)

        ax.text(-0.01, 0.5, ylabel, fontsize=8, rotation=0, va='center', ha='right', transform=ax.transAxes)

        return

    mmin, mmax = np.nanmin(np.nanmin(df.values)), np.nanmax(np.nanmax(df.values))

    """
    if majorMetric == 'correlation':
        missingFillValue = 1.0
    elif majorMetric == 'cosine':
        missingFillValue = 1.0
    elif majorMetric == 'euclidean':
        missingFillValue = mmax
    else:
        missingFillValue = mmax
    """
    missingFillValue = 0
    
    print('Filing missing values with:', missingFillValue, flush=True)
    M = pdist(df.fillna(missingFillValue).values.T, metric=metric)

    fig = plt.figure(figsize=(8, 12))
    dataArgs = addDendro(fig, df.columns, M, [0.1, 0.8, 0.75, 0.165])
    addHeatmap(fig, dataArgs, [0.1, 0.1, 0.75, 0.4])

    print("Making Dendro Dist")
    dendro_dist =[]
    
    if False:
        allGenes = dataArgs['allGenes'] 
        print(len(allGenes))
        for i in range(len(allGenes)):
            dendro_dist.append([])
            for j in range(len(allGenes)):
                g1 = allGenes[i]
                g2 = allGenes[j]
                dendro_dist[i].append(abs(i-j))
        dendro_dist = pd.DataFrame(dendro_dist,allGenes,allGenes)

    print("Plotting")
    
    st, delta = 0.52, 0.035
    addBar(fig, dataArgs, 0, [0.1, st + 0*(delta + 0.015), 0.75, delta])
    addBar(fig, dataArgs, 3, [0.1, st + 1*(delta + 0.015), 0.75, delta])
    addBar(fig, dataArgs, 8, [0.1, st + 2*(delta + 0.015), 0.75, delta])
    #addBar(fig, dataArgs, 11, [0.1, st + 3*(delta + 0.015), 0.75, delta])
    addBar(fig, dataArgs, 12, [0.1, st + 3*(delta + 0.015), 0.75, delta])
    addBar(fig, dataArgs, 11, [0.1, st + 4*(delta + 0.015), 0.75, delta])
    print("plotting last")
    #addBar(fig, dataArgs, 8, [0.1, st + 4*(delta + 0.015), 0.75, delta])
    #print(dataArgs["clusters"])
    #print(len(dataArgs["clusters"]))
    #print(df.shape)
    #print(dataArgs["allGenes"])
    
    
   
    #print(dataArgs["order"])
    fig.suptitle('Data: %s (%s receptors)' % (suffix, df.shape[1]), fontsize=12)
    fig.savefig(os.path.join(saveDir, '%s dendrogram-heatmap-%s.png' % (suffix, majorMetric)), dpi=600)
    plt.close(fig)

    return
