from decneo.commonFunctions import *
from decneo.analysisPipeline import Analysis

if __name__ == '__main__':

    wdir = 'd:/Projects/A_Endothelial/VS/Endothelial/results/L6B6/' 

    PCN = read(wdir + '../../data/PCN')

    def addExtData(ext, extdata, an, common = True):

        def add(panelname, paneldata):

            an.panels += [panelname]
            extdata.update({panelname: paneldata})

            return

        if common:
            extdata.update({'conservedGenes': pd.read_excel(os.path.join(an.bootstrapDir, 'All', 'comparison.xlsx'), index_col=1, header=0)['Inter-measures.T50_common_count']})

            # Bootstrap frequency ---------------------------------------------------------------------------------
            add('Bootstrap of 3', pd.read_excel(an.workingDir + 'Avg combo3avgs_variant.xlsx', index_col=0, header=0)['Bootstrap'])
            add('Bootstrap of 4', pd.read_excel(an.workingDir + 'Avg combo4avgs_variant.xlsx', index_col=0, header=0)['Bootstrap'])

            return extdata

        # TMM ------------------------------------------------------------------------------------------------
        receptors = pd.read_excel(os.path.join(an.bootstrapDir, 'All', 'per-gene-measures-correlation.xlsx'), index_col=[0,1], header=0).xs('Dm', level='approach').index.values

        dfTMM = pd.read_csv(wdir + '%s_PL_TMM.csv' % ext, index_col=0, header=0)
        seB = dfTMM[dfTMM.columns[0:1]].mean(axis=1).replace(0.,np.nan).dropna()
        seL = dfTMM[dfTMM.columns[2:3]].mean(axis=1).replace(0.,np.nan).dropna()
        seLbyB = seL.reindex(seB.index).fillna(0.)/seB
        log2seLbyB = np.log2(seLbyB).replace(np.inf, 0.).replace(-np.inf, 0.).replace(np.nan, 0.)

        add('%s B6 expr.' % ext, seB)
        add('%s L6 expr.' % ext, seL)
        add('%s L6 1st ne. avg. expr.' % ext, pd.Series({receptor:seL.loc[seL.index.intersection(list(PCN.neighbors(receptor)))].mean() for receptor in receptors if receptor in PCN.nodes()}))
        add('%s B6 1st ne. avg. expr.' % ext, pd.Series({receptor:seB.loc[seB.index.intersection(list(PCN.neighbors(receptor)))].mean() for receptor in receptors if receptor in PCN.nodes()}))
        add('%s log2(L6/B6)' % ext, log2seLbyB)

        # DEG ------------------------------------------------------------------------------------------------
        dfDEG = pd.read_excel(wdir + 'df_%s_DEG.xlsx' % ext, index_col=0, header=0)
        
        add('%s Is DEG?' % ext, dfDEG['up'])

        # Network enrichment DEG ------------------------------------------------------------------------------
        add('%s Net. enrich. DEG' % ext, -np.log(binomialEnrichmentProbability('PCN.txt', enriched_genes=dfDEG['padj'].loc[dfDEG['up']==1].index.values, target_genes=receptors, PCNpath=an.PCNpath)['Binomial_Prob']))

        return extdata
        
    kwargs = dict(genesOfInterest=receptorsListHugo_2555, knownRegulators=gEC22 + ['LIFR', 'IL6ST'], 
                  panels = ['fraction', 'binomial', 'markers', 'top50', 'combo3avgs', 'combo4avgs'])

    if False:
        print('Voigt choroid dataset')
        an = Analysis(workingDir=wdir + 'choroid Voigt remapped/', otherCaseDir=wdir + 'fromPanglaoDBmouseAllbyDCS/', **kwargs)
    else:
        print('PanglaoDB mouse dataset')
        an = Analysis(workingDir=wdir + 'PangaloDB_byDCS_mouse/', otherCaseDir=wdir + 'fromPanglaoDBhumanAllbyDCS/', **kwargs)

    customData = externalPanelsData.copy()
    customData = addExtData(None, customData, an)
    customData = addExtData('BCE', customData, an, common=False)
    customData = addExtData('BAE', customData, an, common=False)

    an.analyzeCase(None, dpi=600, suffix='All', saveDir=os.path.join(an.bootstrapDir, 'All/'), toggleCalculateMajorMetric=False, toggleCalculateMeasures=False, externalPanelsData=customData, includeClusterNumber=False, togglePublicationFigure=True, toggleExportFigureData=False, toggleIncludeHeatmap=True, toggleAdjustText=True)