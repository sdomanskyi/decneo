from decneo.commonFunctions import *
import DigitalCellSorter
from decneo.analysisPipeline import Analysis, process

if platform.system()=="Windows":
    dataDir = 'results/new_validation/' 
    wdir = 'results/new_validation/' 
else:
    dataDir = '/mnt/research/piermarolab/new_validation/'
    wdir = '/mnt/research/piermarolab/Sergii/secondValidation/'

if __name__ == '__main__':

    print('Current working directory:', os.getcwd(), flush=True)

    # Prepare fibroblasts
    if False:
        datasetName = 'E-MTAB-7149'
        DCS = DigitalCellSorter.DigitalCellSorter(saveDir=wdir + datasetName, dataName=datasetName, species='Mouse')

        fileName = datasetName + '.h5'
        filePath = dataDir + datasetName + '/' + fileName
        df = pd.read_hdf(filePath, key=KeysOfStore(filePath)[0]).fillna(0.)
        tempIndex = df.index.values.copy()
        tempIndex[tempIndex==tempIndex] = np.array(DCS.gnc.Convert(tempIndex[tempIndex==tempIndex].tolist(), 'ensembl', 'hugo', returnUnknownString=False))
        df.index = pd.Series(tempIndex).replace(Mouse_to_Human_HUGO_conversion).values
        df.columns = df.columns.str.split('-', n=1, expand=True)
        df.columns.names = ['batch', 'cell']
        print(df)

        cells = DCS.getCells(clusterIndex='2.0')
        print(len(cells))

        df = df[cells]
        df = df.loc[df.sum(axis=1) > 0.]
        print(df)

        df.to_hdf(wdir + 'E-MTAB-7149_3536-Fibroblasts.h5', key='df', **phdf)

    # Annotate with DCS and extract EC/nonEC cells
    if False:
        newDatasetsInfo = pd.read_excel(wdir + 'Study_Summary.xlsx', index_col=0, header=0)
        #newDatasetsInfo = newDatasetsInfo.loc[newDatasetsInfo['Ready to process']==1]
        newDatasetsInfo = newDatasetsInfo.loc[newDatasetsInfo['Processed']==1]

        print(newDatasetsInfo)

        dfsh = []
        dfsm = []
        dfshNon = []
        dfsmNon = []
        for datasetName in newDatasetsInfo.index.values:
            species = newDatasetsInfo.loc[datasetName]['Species']
            print('\n', datasetName, '\t', species, flush=True)

            try:
                DCS = DigitalCellSorter.DigitalCellSorter(saveDir=wdir + datasetName, dataName=datasetName, species=species, geneNamesType='alias', doBatchCorrection=False, precutQC=True, nClusters=12, mitochondrialGenesCutoffQC=1.75, geneListFileName='/mnt/home/domansk6/Projects/Endothelial/dev/geneLists/fibro_endo_epi_v2_Human.xlsx', updateConversionDictFile=False, layout='TSNE', availableCPUsCount=10, verbose=0)

                if False:
                    batchSeparator = newDatasetsInfo.loc[datasetName]['Batch separator']
                    fileName = [fileName for fileName in os.listdir(dataDir + datasetName) if datasetName in fileName and fileName[-3:]=='.h5'][0]
                    filePath = dataDir + datasetName + '/' + fileName
                    key = KeysOfStore(filePath)[0]
                    df = pd.read_hdf(filePath, key=key).fillna(0.)

                    if datasetName == 'GSE155788':
                        df.index = df.index.str.split('\t', n=1, expand=True).get_level_values(0)

                    fromType = 'ensembl'
               
                    if datasetName == 'GSE159354' or datasetName == 'GSE156410':
                        fromType = 'alias'   

                    tempIndex = df.index.values.copy()
                    tempIndex[tempIndex==tempIndex] = np.array(DCS.gnc.Convert(tempIndex[tempIndex==tempIndex].tolist(), fromType, 'hugo', returnUnknownString=False))
                    df.index = tempIndex

                    df.index = pd.Series(df.index.values).replace(Mouse_to_Human_HUGO_conversion).values

                    if not batchSeparator=='None':
                        df.columns = df.columns.str.split(batchSeparator, n=1, expand=True)
                        df.columns.names = ['batch', 'cell']

                    if datasetName == 'E-MTAB-7427':
                        df.columns.names = ['cell']
                        df = pd.concat([df], axis=1, keys=['E-MTAB-7427'], names=['batch'], sort=False).reorder_levels(['batch', 'cell'], axis=1)
                        print(df)

                    df = df.groupby(level=0, axis=0).sum()
                    df.to_hdf(wdir + datasetName + '.h5', key='df', **phdf)

                    if True:
                        for gene in ['ENSMUSG00000062960','ENSMUSG00000020717','Kdr','Pecam1','ENSG00000128052','ENSG00000261371','KDR','PECAM1']:
                            try:
                                se = df.loc[gene].fillna(0.)
                                fr = (se > 0.).sum() / len(se)
                                print(gene, '\t', len(se), '\t', (se > 0.).sum(), '\t', np.round(fr, 2), '\t', flush=True)
                            except Exception as exception:
                                pass

                    if datasetName != 'E-MTAB-7149':
                        dff = pd.read_hdf(wdir + 'E-MTAB-7149_3536-Fibroblasts.h5', key='df')
                        dff = dff.groupby(level=0, axis=0).sum()
                        df = pd.concat([df, dff], axis=1, sort=False).fillna(0.)
                        df = df.T.loc[~df.T.index.duplicated(keep='first')].T

                        DCS.excludedFromQC = df.columns[np.isin(df.columns.get_level_values('batch').values, np.array(['ERR2737398', 'ERR2737399', 'ERR2737400']))]
                        print(len(DCS.excludedFromQC), flush=True)

                    DCS.prepare(df)

                    del df

                    print(DCS.df_expr, flush=True)

                    DCS.calculateQCmeasures()
                    DCS.normalize(median=DCS.medianScaleFactor)
                    DCS.project()
                    DCS.qualityControl()
                    DCS.cluster()
                    DCS.recordExpressionData()
                    DCS.annotate()
                    DCS.makeAnnotationResultsMatrixPlot()
                    DCS.makeMarkerExpressionPlot()
                    DCS.makeStackedBarplot()
                    DCS.makeProjectionPlotAnnotated()

                if True:
                    EC = DCS.getCells('Endothelial cells')
                    
                    if datasetName != 'E-MTAB-7149':
                        EC = EC.difference(pd.read_hdf(wdir + 'E-MTAB-7149_3536-Fibroblasts.h5', key='df').columns)

                    df_temp = DCS.getExprOfCells(EC).droplevel('cluster', axis=1)
                    print(df_temp.shape, flush=True)

                    df_temp.columns.set_levels(datasetName + ' ' + df_temp.columns.levels[0], level=0, inplace=True)

                    if species=='Mouse':
                        dfsm.append(df_temp)
                    elif species=='Human':
                        dfsh.append(df_temp)
                    else:
                        print(species, flush=True)

                if True:
                    nonECall = DCS.getCells().index.difference(DCS.getCells('Endothelial cells'))
                    print('All non-EC (with FB): %s' % nonECall.shape, end='\t')

                    if datasetName != 'E-MTAB-7149':
                        nonECall = nonECall.difference(pd.read_hdf(wdir + 'E-MTAB-7149_3536-Fibroblasts.h5', key='df').columns)
                        print('All non-EC: %s' % nonECall.shape, end='\t')

                    if True:
                        dh = dict(zip(*np.unique(nonECall.get_level_values('batch').values, return_counts=True)))

                        temps = []
                        for key in dh.keys():
                            #print(key, '\t', dh[key])
                            temp = nonECall[nonECall.get_level_values('batch').values == key]
                            temp = pd.Series(index=temp, data=0).sample(n=min(1000, len(temp)), replace=False)
                            temps.append(temp)

                        nonEC = pd.concat(temps, axis=0, sort=False).index
                    else:
                        nonEC = pd.Series(index=nonECall, data=0).sample(n=min(1000, len(nonECall)), replace=False).index
                    
                    print('Selected non-EC: %s' % nonEC.shape, end='\n')

                    df_temp = DCS.getExprOfCells(nonEC).droplevel('cluster', axis=1)
                    print('Data shape:', df_temp.shape, flush=True)

                    df_temp.columns.set_levels(datasetName + ' ' + df_temp.columns.levels[0], level=0, inplace=True)

                    if species=='Mouse':
                        dfsmNon.append(df_temp)
                    elif species=='Human':
                        dfshNon.append(df_temp)
                    else:
                        print(species, flush=True)

            except Exception as exception:
                print(exception)

        if True:
            dfsm = pd.concat(dfsm, axis=1, sort=False).fillna(0.)
            print(dfsm, flush=True)
            dfsm.to_hdf(wdir + 'secondValidationEC.h5', key='df_mouse', **phdf)

            dfsh = pd.concat(dfsh, axis=1, sort=False).fillna(0.)
            print(dfsh, flush=True)
            dfsh.to_hdf(wdir + 'secondValidationEC.h5', key='df_human', **phdf)

        if True:
            dfsmNon = pd.concat(dfsmNon, axis=1, sort=False).fillna(0.)
            print(dfsmNon, flush=True)
            dfsmNon.to_hdf(wdir + 'secondValidationNonEC.h5', key='df_mouse', **phdf)

            dfshNon = pd.concat(dfshNon, axis=1, sort=False).fillna(0.)
            print(dfshNon, flush=True)
            dfshNon.to_hdf(wdir + 'secondValidationNonEC.h5', key='df_human', **phdf)

    # DECNEO bootstrap and analysis
    if False:
        anMouse = process(*(None, None), *(None, None),
                        wdir + 'DECNEO analysis/', '/mnt/research/piermarolab/Sergii/results/PanglaoDB_byDCS_human_correlation/', 
                        nCPUs=4 if platform.system()=="Windows" else 20, parallelBootstrap=True,
                        genesOfInterest=receptorsListHugo_2555, knownRegulators=gEC22, perEachOtherCase=True,
                        nBootstrap=100, part1=False, part2=False, part3=False)[0]

        if False:
            dfa = pd.read_hdf(wdir + 'secondValidationEC.h5', key='df_mouse')
            print(dfa, flush=True)

            dfb = pd.read_hdf(wdir + 'secondValidationNonEC.h5', key='df_mouse').reindex(dfa.index).fillna(0.)
            print(dfb, flush=True)

            anMouse.prepareDEG(dfa, dfb)
            anMouse.preparePerBatchCase(exprCutoff=0.05)
        
        anMouse.prepareBootstrapExperiments(parallel=True)
        anMouse.analyzeBootstrapExperiments()
        anMouse.reanalyzeMain()
        anMouse.analyzeCombinationVariant('Avg combo3avgs')
        anMouse.analyzeCombinationVariant('Avg combo4avgs')
        anMouse.reanalyzeMain()
        anMouse.compareTwoCases(wdir + 'DECNEO analysis/bootstrap/All/', '/mnt/research/piermarolab/Sergii/results/PanglaoDB_byDCS_human_correlation/bootstrap/All/', saveName=wdir + 'DECNEO analysis/bootstrap/All/comparison')
        anMouse.reanalyzeMain(togglePublicationFigure=True, markersLabelsRepelForce=1.5, includeClusterNumber=False)

    if True:
        #dfa = pd.read_hdf(wdir + 'secondValidationEC.h5', key='df_mouse')
        #print(dfa, flush=True)
        #dfa.loc['KDR'].to_excel(wdir + 'KDR_EC.xlsx')
        #dfa.loc['FLT1'].to_excel(wdir + 'FLT1_EC.xlsx')
        #dfa.loc['PECAM1'].to_excel(wdir + 'PECAM1_EC.xlsx')

        #dfb = pd.read_hdf(wdir + 'secondValidationNonEC.h5', key='df_mouse').reindex(dfa.index).fillna(0.)
        #print(dfb, flush=True)
        #dfb.loc['KDR'].to_excel(wdir + 'KDR_nonEC.xlsx')
        #dfb.loc['FLT1'].to_excel(wdir + 'FLT1_nonEC.xlsx')
        #dfb.loc['PECAM1'].to_excel(wdir + 'PECAM1_nonEC.xlsx')

        dir = 'd:/Projects/A_Endothelial/VS\Endothelial/results/for meeting 12 28 2020/EC nonEC/'

        def getSe(gene, kind = ''):

            se = pd.read_excel(dir + '%s_%sEC.xlsx' % (gene, kind), index_col=[0,1], header=0).fillna(0.)

            index3 = pd.Series(index=se.index.get_level_values('batch').str.split(' ', expand=True), data=se.index.get_level_values('cell')).to_frame().set_index('cell', append=True).index
            index3.names = ['study', 'sample', 'cell']

            se.index = index3

            return se

        temp = getSe('KDR', kind='').groupby(level=[0, 1]).count()
        temp = temp.loc[temp.values>=10]
        studiesIndex = temp.index
        #getSe('KDR', kind='non').groupby(level=[0, 1]).count().reindex(studiesIndex).to_excel(dir + 'gNon.xlsx', merge_cells=False)

        #print(getSe(gene).groupby(level=[0, 1]).count())

        for gene in ['KDR', 'FLT1', 'PECAM1']:
            se = getSe(gene, kind='non').replace(0., np.nan).groupby(level=[0,1]).count().reindex(studiesIndex).groupby(level=0).sum()
            se.to_excel(dir + 'gs.xlsx', merge_cells=False)
            input('***')
