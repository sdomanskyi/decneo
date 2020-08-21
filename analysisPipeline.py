from common import *

def analyze(df_expr, selGenes, stimulators, inhibitors, majorMetric, toggleCalculateMajorMetric = True, exprCutoff = 0.05, toggleExportFigureData = True, toggleCalculateMeasures = True, suffix = '', saveDir = '', toggleGroupBatches = True, dpi = 300, toggleAdjustText = True, panels = None, figureSize=(8, 22), toggleAdjustFigureHeight=True, noPlot = False, halfWindowSize = 10, printStages = True, externalPanelsData = None, metricsFile = 'metricsFile.h5', toggleIncludeHeatmap = True):

    '''Parameters:
        df_expr: Take one species, one cluster (subset of clusters)
        selGenes: List of receptors
        stimulators: will be highlighted on plots in green
        inhibitors: will be highlighted on plots in red
        metric: major matric used
        saveDir: exerything is exported to this directory, should be unique for each dataset
    '''

    np.random.seed(0)

    def calculateMajorMetricAndGeneStats(df_expr, majorMetric, saveDir, groupBatches, selGenes, exprCutoff):

        '''Calculate cdist of metric (e.g. correlation) (median across batches if groupBatches is True).
        Calculate fraction of cells expressing each gene, and median of non-zero gene expression (per batch)'''

        print('Received expression data of shape:', df_expr.shape, flush=True)
        np.savetxt(os.path.join(saveDir, 'size.txt'), np.array(df_expr.shape), fmt='%i')

        # For each batch calculate gene expression distance metric
        print('Calculating distance metric', flush=True)
        df_measure = get_df_distance(df_expr, metric=majorMetric, genes=selGenes, analyzeBy='batch', minSize=10, groupBatches=groupBatches, cutoff=exprCutoff)

        print('Recording major metric (shape: %s, %s) to h5' % df_measure.shape, flush=True)
        df_measure.to_hdf(os.path.join(saveDir, metricsFile), key=majorMetric, mode='a', complevel=4, complib='zlib')

        # For each batch calculate fraction of cells expressing each gene
        df_fraction = df_expr.replace(0, np.nan).replace(0., np.nan).groupby(axis=1, level='batch').agg('count') /\
            df_expr.fillna(0.).groupby(axis=1, level='batch').agg('count')

        print('Recording fractions (shape: %s, %s) to h5' % df_fraction.shape, flush=True)
        df_fraction.to_hdf(os.path.join(saveDir, 'perGeneStats.h5'), key='df_fraction', mode='a', complevel=4, complib='zlib')

        # For each batch calculate median of non-zero values of each gene expression
        df_median_expr = df_expr.replace(0, np.nan).replace(0., np.nan).groupby(axis=1, level='batch').agg(np.nanmedian)

        print('Recording median expression (shape: %s, %s) to h5' % df_fraction.shape, flush=True)
        df_median_expr.to_hdf(os.path.join(saveDir, 'perGeneStats.h5'), key='df_expression', mode='a', complevel=4, complib='zlib')

        # For each batch calculate median of non-zero values of each gene expression
        se_count = df_expr.fillna(0.).groupby(axis=1, level='batch').agg('count').iloc[0]

        print('Recording per batch counts to h5', flush=True)
        se_count.to_hdf(os.path.join(saveDir, 'perGeneStats.h5'), key='se_count', mode='a', complevel=4, complib='zlib')

        return

    def makeCombinationPlot(df, metric = 'euclidean', linkageMethod = 'ward', n_clusters = 10, adjustText = toggleAdjustText):

        '''The metric is to build dendrogram and identify clusters in it.
        Metric has to be of type "Euclidean" to use linkage method "Ward".
        With any other metric (e.g. correlation distance) use linkage method "average" etc.

        metric 'euclidean_missing' used commonly-non-missing points only
        '''
             
        nonlocal panels, figureSize

        if panels is None:
            panels = [
                    #'combo3', 
                    'combo3avgs', 
                    
                    #'variability_3mean',                    
                    #'variability_3std',                    
                    #'variability_3cov',

                    #'combo4', 
                    #'combo4avgs',  
                    
                    #'variability_4mean',                    
                    #'variability_4std',                    
                    #'variability_4cov',   
                    
                    'fraction', 
                    'binomial', 
                    'top50', 
                    #'markers', 

                    #'PubMedHits', 
                    #'gAbove50_PanglaoMouse', 
                    #'gAbove50_PanglaoHuman', 
                    #'GOpositive', 
                    #'GOnegative', 

                    #'markerstop50', 
                    'expression', 
                    'closeness',
                    'age', 
                    'rate', 
                    ]

        def addDendro(fig, dataGenes, M, coords, metric=metric, linewidth=0.25, adjustText = adjustText, fontsize = 5):

            genesSubset = list(stimulators) + list(inhibitors)

            ax = fig.add_axes(coords, frame_on=False)

            Z = hierarchy.linkage(np.nan_to_num(M, nan=max(M)), method='ward', optimal_ordering=True)

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
                        texts.append(ax.text(xpos, -2., xtext, fontsize=fontsize, rotation=90, va='top', ha='center', color=color))
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
                    ltext = ax.text(position, vposition, '#%s' % cluster, fontsize=fontsize, color='white', va='center', ha='center')
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

        def addHeatmap(fig, dataArgs, coords, adjustText = adjustText, fontsize = 5):

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
                        texts.append(ax.text(xpos, 1.01*ax.get_ylim()[0], xtext, fontsize=fontsize, rotation=90, va='top', ha='center', color=color))
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
                        texts.append(ax.text(-0.01*ax.get_xlim()[1], ypos, xtext, fontsize=fontsize, va='center', ha='right', color=color))
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
                    ltext = ax.text(position, position, '#%s' % cluster, fontsize=fontsize, color='white', va='center', ha='center')
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
                ax = fig.add_axes([0.85, 0.05, 0.025, 0.1], frame_on=False)

                ax.set_xticks([])
                ax.set_xticklabels([])
                ax.set_yticks([])
                ax.set_yticklabels([])

                clb = fig.colorbar(im, ax=ax, fraction=0.4, label='Eucl. dist. of gene expr.\n %s dist.' % majorMetric)
                clb.ax.tick_params(labelsize=fontsize)

            return

        def addBar(fig, dataArgs, panel, coords, halfWindowSize = halfWindowSize):

            nonlocal panelsData, panelsDataNames, externalPanelsData

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

            if panel == 'fraction':
                ylabel='Fraction'
                try:
                    data =pd.read_hdf(os.path.join(saveDir, 'perGeneStats.h5'), key='df_fraction').reindex(allGenes)
                    if type(data) is pd.Series:
                        data = data.values
                    else:
                        data = data.mean(axis=1).values
                except:
                    data = np.zeros(len(allGenes))

            elif panel == 'expression':
                ylabel='Expression'
                try:
                    data = pd.read_hdf(os.path.join(saveDir, 'perGeneStats.h5'), key='df_expression').reindex(allGenes)
                    if type(data) is pd.Series:
                        data = data.values
                    else:
                        data = data.mean(axis=1).values
                except:
                    data = np.zeros(len(allGenes))

            elif panel == 'closeness':
                ylabel='Closeness'
                try:
                    data = -np.append([0], np.diagonal(M, offset=1))
                    data -= min(data)
                except:
                    data = np.zeros(len(allGenes))

            elif panel == 'rate':
                ylabel='Rate'
                try:
                    data = externalPanelsData['Evolutionary rate'].reindex(allGenes).values              
                except:
                    data = np.zeros(len(allGenes))

            elif panel == 'age':
                ylabel='Age'
                try:
                    data = externalPanelsData['Evolutionary age'].reindex(allGenes).values 
                except:
                    data = np.zeros(len(allGenes))

            elif panel == 'markerstop50':
                ylabel='Markers\nTop50\noverlap'
                try:
                    data = externalPanelsData['conservedMarkers']
                    data = data.loc[~data.index.duplicated(keep='first')].reindex(allGenes).values
                except:
                    data = np.zeros(len(allGenes))

            elif panel == 'PubMedHits':
                ylabel='log(PubMedHits)'
                try:
                    data = np.log(pd.Series(externalPanelsData['pubMed angiogenesis hits'])).reindex(allGenes).values
                except:
                    data = np.zeros(len(allGenes))

            elif panel == 'gAbove50_PanglaoMouse':
                ylabel='gAbove50\nPanglaoMouse'
                try:
                    data = np.zeros(len(allGenes))
                    data[np.where(np.isin(allGenes, externalPanelsData['gAbove50_PanglaoMouse']))[0]] = 1.
                except:
                    data = np.zeros(len(allGenes))

            elif panel == 'GOpositive':
                ylabel='GO positive'
                try:
                    data = np.zeros(len(allGenes))
                    data[np.where(np.isin(allGenes, externalPanelsData['GOpositive']))[0]] = 1.
                except:
                    data = np.zeros(len(allGenes))

            elif panel == 'GOnegative':
                ylabel='GO negative'
                try:
                    data = np.zeros(len(allGenes))
                    data[np.where(np.isin(allGenes, externalPanelsData['GOnegative']))[0]] = 1.
                except:
                    data = np.zeros(len(allGenes))

            elif panel == 'gAbove50_PanglaoHuman':
                ylabel='gAbove50\nPanglaoHuman'
                try:
                    data = np.zeros(len(allGenes))
                    data[np.where(np.isin(allGenes, externalPanelsData['gAbove50_PanglaoHuman']))[0]] = 1.
                except:
                    data = np.zeros(len(allGenes))

            elif panel == 'markers':
                ylabel='Markers'
                try:
                    data = np.zeros(len(allGenes))
                    data[locations] = 1.
                except:
                    data = np.zeros(len(allGenes))

            elif panel == 'top50':
                ylabel='Top50\noverlap'
                try:
                    data = externalPanelsData['conservedGenes']
                    data = data.loc[~data.index.duplicated(keep='first')].reindex(allGenes).values
                except:
                    data = np.zeros(len(allGenes))

            elif panel == 'binomial':
                ylabel='Binomial\n-log(pvalue)'
                try:
                    data = pd.read_hdf(os.path.join(saveDir, 'perGeneStats.h5'), key='df_ranks')
                    if not type(data) is pd.Series:
                        data = data.median(axis=1)

                    data = data[data > 0].sort_values()[:1000]
                    diffExpressedGenes = data.index
                    data = -np.log(nx_binom('data/PCN.txt', enriched_genes=diffExpressedGenes, target_genes=receptorsListHugo_2555)['Binomial_Prob'].reindex(allGenes).values)
                except:
                    data = np.zeros(len(allGenes))

                    try:
                        data = externalPanelsData['diffExpressedGenes']
                        if not type(data) is pd.Series:
                            data = data.median(axis=1)

                        data = data[data > 0].sort_values()[:1000]
                        diffExpressedGenes = data.index

                        data = -np.log(nx_binom('data/PCN.txt', enriched_genes=diffExpressedGenes, target_genes=receptorsListHugo_2555)['Binomial_Prob'].reindex(allGenes).values)
                    except:
                        pass

            elif panel == 'variability_3mean':
                ylabel='variability\n3 mean'
                try:
                    data = externalPanelsData['variability_3']['mean'].reindex(allGenes).values              
                except:
                    data = np.zeros(len(allGenes))

            elif panel == 'variability_3std':
                ylabel='variability\n3 std'
                try:
                    data = externalPanelsData['variability_3']['std'].reindex(allGenes).values              
                except:
                    data = np.zeros(len(allGenes))

            elif panel == 'variability_3cov':
                ylabel='variability\n3 cov'
                try:
                    data = externalPanelsData['variability_3']['cov'].reindex(allGenes).values              
                except:
                    data = np.zeros(len(allGenes))
                    
            elif panel == 'variability_4mean':
                ylabel='variability\n4 mean'
                try:
                    data = externalPanelsData['variability_4']['mean'].reindex(allGenes).values              
                except:
                    data = np.zeros(len(allGenes))

            elif panel == 'variability_4std':
                ylabel='variability\n4 std'
                try:
                    data = externalPanelsData['variability_4']['std'].reindex(allGenes).values              
                except:
                    data = np.zeros(len(allGenes))

            elif panel == 'variability_4cov':
                ylabel='variability\n4 cov'
                try:
                    data = externalPanelsData['variability_4']['cov'].reindex(allGenes).values              
                except:
                    data = np.zeros(len(allGenes))

            elif panel == 'combo3':
                ylabel='Combination\n3'
                try:
                    data = np.zeros(len(allGenes))

                    def func(s):
                       
                        w = np.nansum(panelsData[panelsDataNames[s]])
                        if w != w or w == 0.:
                            w = 1.

                        result = np.nan_to_num(panelsData[panelsDataNames[s]]) / w

                        return result

                    data += func('fraction')
                    data += func('binomial')
                    data += func('top50')
                except:
                    data = np.zeros(len(allGenes))

            elif panel == 'combo3avgs':
                ylabel='Combination\n3 avgs'
                try:
                    data = np.zeros(len(allGenes))

                    def func(s):
                       
                        w = np.nansum(panelsData[panelsDataNames[s]])
                        if w != w or w == 0.:
                            w = 1.

                        result = np.nan_to_num(panelsData[panelsDataNames[s]]) / w

                        return movingAverageCenteredLooped(result, halfWindowSize)

                    data += func('fraction')
                    data += func('binomial')
                    data += func('top50')
                except:
                    data = np.zeros(len(allGenes))
                    
            elif panel == 'combo4':
                ylabel='Combination\n4'
                try:
                    data = np.zeros(len(allGenes))

                    def func(s):
                       
                        w = np.nansum(panelsData[panelsDataNames[s]])
                        if w != w or w == 0.:
                            w = 1.

                        result = np.nan_to_num(panelsData[panelsDataNames[s]]) / w

                        return result

                    data += func('fraction')
                    data += func('binomial')
                    data += func('top50')
                    data += func('markers')
                except:
                    data = np.zeros(len(allGenes))

            elif panel == 'combo4avgs':
                ylabel='Combination\n4 avgs'
                try:
                    data = np.zeros(len(allGenes))

                    def func(s):
                       
                        w = np.nansum(panelsData[panelsDataNames[s]])
                        if w != w or w == 0.:
                            w = 1.

                        result = np.nan_to_num(panelsData[panelsDataNames[s]]) / w

                        return movingAverageCenteredLooped(result, halfWindowSize)

                    data += func('fraction')
                    data += func('binomial')
                    data += func('top50')
                    data += func('markers')
                except:
                    data = np.zeros(len(allGenes))

            ax.bar(range(len(clusters)), data, width=ax.get_xlim()[1]/len(clusters), color=tickLabelsColors)

            data_avg = movingAverageCenteredLooped(np.nan_to_num(data), halfWindowSize)
            ax.plot(range(len(clusters)), data_avg, linewidth=1.0, color='coral', alpha=1.0)

            ax.text(0.999, 0.95, 'window = %s' % (2*halfWindowSize + 1), c='coral', ha='right', va='top', transform=ax.transAxes, fontsize=3)

            nandata = np.isnan(data)
            if np.sum(nandata) > 0:
                wh = np.where(nandata)
                ax.bar(wh[0], [np.nanmax(data)]*wh[0].shape[0], width=ax.get_xlim()[1]/len(clusters), color='lightgrey', alpha=0.3)

            ax.set_xticks([])
            ax.set_xticklabels([])

            yticks = np.round(ax.get_ylim(), 3)
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticks)

            ax.tick_params(axis='y', labelsize=3.5, width=0.75, length=3)

            if True:
                ylim = ax.get_ylim()
                for i in range(1, len(np.unique(clusters))):
                    ax.plot([clusterBoundaries[i] - 0.5]*2, [ylim[0], ylim[1]], '--', lw=0.5, color='k', clip_on=False)

            ax.text(-0.01, 0.5, ylabel, fontsize=4, rotation=0, va='center', ha='right', transform=ax.transAxes)

            return ylabel, data, data_avg

        if metric == 'euclidean_missing':
            metric = metric_euclidean_missing

        mmin, mmax = np.nanmin(np.nanmin(df.values)), np.nanmax(np.nanmax(df.values))

        if majorMetric == 'correlation':
            missingFillValue = 1.0
        elif majorMetric == 'cosine':
            missingFillValue = 1.0
        elif majorMetric == 'euclidean':
            missingFillValue = mmax
        else:
            missingFillValue = mmax

        if printStages:
            print('\tCalculating %s metric of %s . . .' % (metric, majorMetric), end='\t', flush=True)
        M = pdist(df.fillna(missingFillValue).values.T, metric=metric)
        if printStages:
            print('Done', flush=True)

        nPanels = len(panels)
        dendroHeight = 0.10
        panelHeight = 0.022
        detla = 0.01

        vsum = nPanels*(panelHeight + detla) + dendroHeight + 0.8

        if toggleAdjustFigureHeight:
            figureSize = (figureSize[0],  figureSize[0] * (9./8.) * vsum)

        factor = 1 / vsum

        topBorder = 1. - 0.03*factor
        bottomBorder = 0.05*factor

        dendroHeight *= factor
        panelHeight *= factor
        detla *= factor

        fig = plt.figure(figsize=figureSize)

        dataArgs = addDendro(fig, df.columns, M, [0.1, topBorder-dendroHeight, 0.75, dendroHeight], metric=metric)

        heatmapHeight = (topBorder - bottomBorder - dendroHeight) - nPanels * (panelHeight + detla) - 0.05*factor
        if nPanels > 0:
            if printStages:
                print('\tPlotting bar panels . . .', end='\t', flush=True)
            panelsData = dict()
            panelsDataNames = dict()
            for ipanel, panel in enumerate(reversed(panels)):
                panelName, data, data_avg = addBar(fig, dataArgs, panel, [0.1, 0.015*factor + bottomBorder + heatmapHeight + ipanel*(panelHeight + detla), 0.75, panelHeight])

                wname = panelName.replace('\n', ' ')

                panelsData.update({wname: data})
                panelsDataNames.update({panel: wname})

                panelsData.update({'Avg ' + wname: data_avg})
                panelsDataNames.update({'Avg ' + panel: 'Avg ' + wname})

            dataArgs.update({'panelsData': panelsData})
            if printStages:
                print('Done', flush=True)

        if not noPlot:
            if toggleIncludeHeatmap:
                if printStages:
                    print('\tPlotting heatmap . . .', end='\t', flush=True)
                addHeatmap(fig, dataArgs, [0.1, bottomBorder, 0.75, heatmapHeight])
                if printStages:
                    print('Done', flush=True)

            number_of_cells = np.loadtxt(os.path.join(saveDir, 'size.txt'), dtype=int)[1]
            fig.suptitle('Data: %s\n(%s cells, %s receptors)' % (suffix, number_of_cells, df.shape[1]), fontsize=8, backgroundcolor='white')

            if printStages:
                print('\tSaving image . . .', end='\t', flush=True)
            fig.savefig(os.path.join(saveDir, '%s dendrogram-heatmap-%s.png' % (suffix, majorMetric)), dpi=dpi)
            if printStages:
                print('Done', flush=True)

        plt.close(fig)

        return dataArgs

    def exportFigureData(dataArgs, saveXLSX = True, saveHDF = True):

        M = dataArgs['M'] 
        allGenes = dataArgs['allGenes'] 
        clusters =  dataArgs['clusters'] 
        selGenes = dataArgs['genes'] 
        selGenesLocations = dataArgs['locations'] 

        df_C = pd.DataFrame(index=allGenes)
        df_C['cluster'] = clusters

        try:
            for panel, panelData in dataArgs['panelsData'].items():
                df_C[panel] = panelData
        except:
            pass

        df_C['stimulator'] = np.where(np.isin(allGenes, stimulators), True, np.nan)
        df_C['inhibitor'] = np.where(np.isin(allGenes, inhibitors), True, np.nan)

        df_M = pd.DataFrame(data=M, index=allGenes, columns=allGenes)

        if saveXLSX:
            writer = pd.ExcelWriter(os.path.join(saveDir, 'dendrogram-heatmap-%s-data.xlsx' % majorMetric))
            df_C.to_excel(writer, 'Cluster index')
            df_M.to_excel(writer, 'Expression distance measure')
            writer.save()

        if saveHDF:
            df_C.to_hdf(os.path.join(saveDir, 'dendrogram-heatmap-%s-data.h5' % majorMetric), key='df_C', mode='a', complevel=4, complib='zlib')
            df_M.to_hdf(os.path.join(saveDir, 'dendrogram-heatmap-%s-data.h5' % majorMetric), key='df_M', mode='a', complevel=4, complib='zlib')

        return

    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    selGenes = np.unique(list(selGenes) + list(stimulators) + list(inhibitors))

    if toggleCalculateMajorMetric and (not df_expr is None):
        calculateMajorMetricAndGeneStats(df_expr, majorMetric, saveDir, toggleGroupBatches, selGenes, exprCutoff)

    # Load and prepare df_measure
    if True:
        try:
            df_measure = pd.read_hdf(os.path.join(saveDir, metricsFile), key=majorMetric)

            if not toggleGroupBatches:
                df_measure = pd.Series(data=np.nanmedian(df_measure.values, axis=1), index=df_measure.index).unstack(0)
        except Exception as exception:
            print(exception, flush=True)

            return

        df_measure = df_measure[df_measure.columns.intersection(np.unique(selGenes))]
        #print(df_measure.shape, flush=True)
    
        for gene in df_measure.columns:
            df_measure.loc[gene, gene] = np.nan

        df_measure = df_measure[df_measure.columns[df_measure.count(axis=0) > 0]]

    # Make dendrogram with heatmap, label clusters, export all(!) figure data and clusters 
    if True:
        dataArgs = makeCombinationPlot(df_measure)
        
        if toggleExportFigureData:
            if printStages:
                print('\tExporting data . . .', end='\t', flush=True)
            exportFigureData(dataArgs)
            if printStages:
                print('Done', flush=True)

    # For each Approach [1-4] calculate (1) AUC of EC23, (2) T50, and (3) EC23 of T50
    if toggleCalculateMeasures:
        if printStages:
            print('\tCalculating measures . . .', end='\t', flush=True)
        df_Em = pd.read_hdf(os.path.join(saveDir, 'dendrogram-heatmap-%s-data.h5' % majorMetric), key='df_M')
        orderedGenes = df_Em.columns
        data_Cm = squareform(pdist(df_measure[orderedGenes].fillna(0.).values.T, metric='correlation'))
        df_Cm = pd.DataFrame(index=orderedGenes, columns=orderedGenes, data=data_Cm)
        df_m = df_measure.loc[orderedGenes, orderedGenes]

        def getForGene(approach, gene, orderedGenes, selGenes23, top = 50):

            if approach == 'm':
                scores = df_m.loc[gene]
            elif approach =='Em':
                scores = df_Em.loc[gene]
            elif approach =='Cm':
                scores = df_Cm.loc[gene]
            elif approach =='Dm':
                scores = np.abs(np.arange(len(orderedGenes)) - np.where(orderedGenes==gene)[0])

            se = pd.Series(index=orderedGenes, data=scores).dropna().sort_values(ascending=True)

            AUC = roc_auc_score(~np.isin(se.index.values, selGenes23), se.values) if len(se) >= 10 else np.nan
            T50 = se.index.values[:top]
            subT50 = np.intersect1d(T50, selGenes23)

            return approach, gene, AUC, len(subT50), cleanListString(T50), cleanListString(subT50)

        list_gEC23 = np.unique(list(stimulators) + list(inhibitors))

        temps = []
        for approach in ['m', 'Cm', 'Em', 'Dm']:
            temps.extend([getForGene(approach, gene, orderedGenes, list_gEC23) for gene in orderedGenes])

        temps = np.array(temps).T
        df = pd.DataFrame(index=pd.MultiIndex.from_arrays([temps[0], temps[1]], names=['approach', 'gene']),
                     data=temps[2:].T, columns=['AUC', 'numEC23T50', 'T50', 'EC23T50'])
        #print(df)

        df.to_excel(os.path.join(saveDir, 'per-gene-measures-%s.xlsx' % majorMetric), merge_cells=False)
        df.to_hdf(os.path.join(saveDir, 'per-gene-measures-%s.h5' % majorMetric), key='df', mode='a', complevel=4, complib='zlib')
        if printStages:
            print('Done', flush=True)

    return

def compareTwoCases(saveDir1, saveDir2, name1 = 'N1', name2='N2', majorMetric = 'correlation', saveName = 'saveName'):

    #name1 = saveDir1.split('/')[-2]
    #name2 = saveDir2.split('/')[-2]

    df1 = pd.read_hdf(os.path.join(saveDir1, 'per-gene-measures-%s.h5' % majorMetric), key='df')
    df2 = pd.read_hdf(os.path.join(saveDir2, 'per-gene-measures-%s.h5' % majorMetric), key='df')

    n23_1 = len(np.intersect1d(np.unique(df1.index.get_level_values('gene').values), gEC23))
    n23_2 = len(np.intersect1d(np.unique(df2.index.get_level_values('gene').values), gEC23))

    #print('Comparing two cases', n23_1, n23_2)

    commonIndex = df1.index.intersection(df2.index)
    df1 = df1.loc[commonIndex]
    df2 = df2.loc[commonIndex]

    df_T50 = pd.concat([df1['T50'].str.replace(' ','').str.split(','), df2['T50'].str.replace(' ','').str.split(',')], keys=[name1, name2], axis=1, sort=False)
    df_EC23T50 = pd.concat([df1['EC23T50'].str.replace(' ','').str.split(','), df2['EC23T50'].str.replace(' ','').str.split(',')], keys=[name1, name2], axis=1, sort=False)
    df_AUC = pd.concat([df1['AUC'].astype(float), df2['AUC'].astype(float)], keys=[name1, name2], axis=1, sort=False)
    df_EC23T50_count = df_EC23T50.applymap(len)

    df_EC23T50_common = df_EC23T50.apply(lambda s: np.intersect1d(s[0], s[1]), axis=1) 
    df_EC23T50_common_count = df_EC23T50.apply(lambda s: len(np.intersect1d(s[0], s[1])), axis=1) 
    df_T50_common = df_T50.apply(lambda s: np.intersect1d(s[0], s[1]), axis=1) 
    df_T50_common_count = df_T50.apply(lambda s: len(np.intersect1d(s[0], s[1])), axis=1) 
    df_AUC_avg = df_AUC.apply(np.mean, axis=1) 

    df_res = pd.concat([df_EC23T50_common.apply(cleanListString), 
                        df_EC23T50_common_count, 
                        df_T50_common.apply(cleanListString), 
                        df_T50_common_count, 
                        #df_AUC_avg, 
                        df1['AUC'].astype(float), 
                        df2['AUC'].astype(float),
                        df1['EC23T50'].str.split(',').apply(len), 
                        df2['EC23T50'].str.split(',').apply(len), 
                        df1['EC23T50'], 
                        df2['EC23T50']], 
                       keys=[('Inter-measures', 'EC23T50_common'), 
                             ('Inter-measures', 'EC23T50_common_count'), 
                             ('Inter-measures', 'T50_common'), 
                             ('Inter-measures', 'T50_common_count'), 
                             #('Inter-measures', 'AUC_avg'), 
                             ('Intra-measures', 'AUC ' + name1), 
                             ('Intra-measures', 'AUC ' + name2),
                             ('Intra-measures', 'EC23 count ' + name1 + ' %s' % n23_1), 
                             ('Intra-measures', 'EC23 count ' + name2 + ' %s' % n23_2), 
                             ('Intra-measures', 'EC23 ' + name1 + ' %s' % n23_1), 
                             ('Intra-measures', 'EC23 ' + name2 + ' %s' % n23_2)],
                       axis=1, sort=False)

    df_res.columns = pd.MultiIndex.from_tuples(df_res.columns)
    
    df_res = df_res.sort_index()
    df_res.to_excel('%s.xlsx' % saveName, merge_cells=False)

    df_res_f = df_res.copy()
    df_res_f = df_res_f.loc[(df_res_f[('Intra-measures', 'EC23 count ' + name1 + ' %s' % n23_1)] >= 5) &
                            (df_res_f[('Intra-measures', 'EC23 count ' + name2 + ' %s' % n23_2)] >= 5) &
                            (df_res_f[('Intra-measures', 'AUC ' + name1)] >= 0.5) &
                            (df_res_f[('Intra-measures', 'AUC ' + name2)] >= 0.5)]
    df_res_f.to_excel('%s_filtered.xlsx' % saveName, merge_cells=False)

    return

def prepareBootstrapExperiments(sourceDir, saveDir, ids = [], majorMetric = 'correlation', df_ranks = None, allData = False):

    print('Reading precalculated data from %s' % sourceDir, flush=True)
    df_measure = pd.read_hdf(os.path.join(sourceDir, 'metricsFile.h5'), key=majorMetric)
    print('df_measure', '%1.1fMb'%(sys.getsizeof(df_measure) / 1024**2), flush=True)

    df_fraction =  pd.read_hdf(os.path.join(sourceDir, 'perGeneStats.h5'), key='df_fraction')
    df_median_expr =  pd.read_hdf(os.path.join(sourceDir, 'perGeneStats.h5'), key='df_expression')
    se_count =  pd.read_hdf(os.path.join(sourceDir, 'perGeneStats.h5'), key='se_count')

    if allData:
        pass
    else:
        for id in ids:
            try:
                saveSubDir = 'Experiment %s' % (id + 1)
                print('\n', saveSubDir, flush=True)
        
                if not os.path.exists(os.path.join(saveDir, saveSubDir)):
                    os.makedirs(os.path.join(saveDir, saveSubDir))

                batches = np.random.choice(df_measure.columns, size=len(df_measure.columns), replace=True)

                print('\tAggregating', flush=True)
                df_measure_temp = pd.Series(data=np.nanmedian(df_measure[batches].values.copy(), axis=1, overwrite_input=True), index=df_measure.index)
        
                print('\tUnstacking', flush=True)
                df_measure_temp = df_measure_temp.unstack(0)

                print('\tRecording', flush=True)
                df_measure_temp.to_hdf(os.path.join(saveDir, saveSubDir, 'metricsFile.h5'), key=majorMetric, mode='a', complevel=4, complib='zlib')

                print('\tPreparing gene stats')
                df_fraction_temp = df_fraction[batches]
                df_fraction_temp.columns = df_fraction_temp.columns + '_' + np.array(range(len(df_fraction_temp.columns))).astype(str)
                df_fraction_temp = df_fraction_temp.mean(axis=1)
                df_fraction_temp.to_hdf(os.path.join(saveDir, saveSubDir, 'perGeneStats.h5'), key='df_fraction', mode='a', complevel=4, complib='zlib')

                df_median_expr_temp = df_median_expr[batches]
                df_median_expr_temp.columns = df_median_expr_temp.columns + '_' + np.array(range(len(df_median_expr_temp.columns))).astype(str)
                df_median_expr_temp = df_median_expr_temp.mean(axis=1)
                df_median_expr_temp.to_hdf(os.path.join(saveDir, saveSubDir, 'perGeneStats.h5'), key='df_expression', mode='a', complevel=4, complib='zlib')

                se_count_temp = se_count[batches]
                se_count_temp.index = se_count_temp.index + '_' + np.array(range(len(se_count_temp.index))).astype(str)
                se_count_temp = se_count_temp.sort_values(ascending=False)
                se_count_temp.to_hdf(os.path.join(saveDir, saveSubDir, 'perGeneStats.h5'), key='se_count', mode='a', complevel=4, complib='zlib')
        
                if not df_ranks is None:
                    df_ranks_temp = df_ranks[batches]
                    df_ranks_temp.columns = df_ranks_temp.columns + '_' + np.array(range(len(df_ranks_temp.columns))).astype(str)
                    df_ranks_temp = df_ranks_temp.median(axis=1).sort_values()
                    df_ranks_temp.to_hdf(os.path.join(saveDir, saveSubDir, 'perGeneStats.h5'), key='df_ranks', mode='a', complevel=4, complib='zlib')

                np.savetxt(os.path.join(saveDir, saveSubDir, 'size.txt'), [df_fraction_temp.shape[0], se_count_temp.sum()], fmt='%i')

            except Exception as exception:
                print(exception)

    return

def runPairOfBootstrapExperiments(args):

    id, saveDir = args
    saveSubDir = 'Experiment %s' % (id + 1)
    print('\n', saveSubDir, flush=True)

    try:
        comparisonName = os.path.join(saveDir, 'Mus musculus', saveSubDir, 'comparison')

        if False:
            for species in ['Homo sapiens', 'Mus musculus']:
                print('Processing:', species)
                analyze(None, receptorsListHugo_2555, gECs, gECi, 'correlation', toggleAdjustText=False, noPlot=True, panels=[],
                        suffix=saveSubDir + ', ' + species, saveDir=os.path.join(saveDir, species, saveSubDir),
                        toggleCalculateMajorMetric=False, toggleExportFigureData=True, toggleCalculateMeasures=True,
                        diffExpressedGenes=None, conservation=None, printStages=False)

            compareTwoCases(os.path.join(saveDir, 'Homo sapiens', saveSubDir, ''), 
                            os.path.join(saveDir, 'Mus musculus', saveSubDir, ''), 
                            name1='human', name2='mouse', saveName=comparisonName)

        if True:
            conservation = {'conservedGenes': pd.read_excel(comparisonName + '.xlsx', index_col=1, header=0)['Inter-measures.T50_common_count'],
                           'conservedMarkers': pd.read_excel(comparisonName + '.xlsx', index_col=1, header=0)['Inter-measures.EC23T50_common_count']}
    
            for species in ['Homo sapiens', 'Mus musculus']:
                print('Re-processing:', species)
                analyze(None, receptorsListHugo_2555, gECs, gECi, 'correlation', toggleAdjustText=False, dpi=300,
                        suffix=saveSubDir + ', ' + species, saveDir=os.path.join(saveDir, species, saveSubDir), 
                        toggleCalculateMajorMetric=False, toggleExportFigureData=True, toggleCalculateMeasures=False,
                        diffExpressedGenes=None, conservation=conservation, printStages=False) 

    except Exception as exception:
        print(exception)

    return
