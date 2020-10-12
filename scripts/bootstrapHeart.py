from commonFunctions import *

# Begin: Part specific for each analyses ----------------------------------------------------------

otherCaseDir = 'results/workflow/PanglaoDB EC all cells w21/Mus musculus/'
workingDir = 'results/Heart_bootstrap_original_annotation/'
nCPUs = 2 if platform.system() == "Windows" else 10

def prepareInputData():

    df = pd.read_hdf('data/GSE109816_heart.h5', key='df')
    origAnno = pd.read_csv('GSE109816.txt', index_col=0, header=0, sep='\t')['CellType']
    origAnno = origAnno.reindex(df.columns.get_level_values('cell')).fillna('FailedQC')
    df.columns = pd.MultiIndex.from_frame(df.columns.to_frame().assign(celltype=origAnno.values))

    import DigitalCellSorter
    Convert = DigitalCellSorter.DigitalCellSorter(matplotlibMode='TkAgg').gnc.Convert
    df.index = Convert(df.index.values.tolist(), 'alias', 'hugo', returnUnknownString=False)
    df = df.loc[~df.index.duplicated(keep='first')]
    df = df.T.loc[~df.T.index.duplicated(keep='first')].T
    df = df.astype(float)
    df = df[df.sum(axis=1) > 0].apply(lambda cell: cell * 10000. / np.sum(cell), axis=0)
    df = np.log2(df.replace(0., np.min(df.values[df.values > 0.])))
    df -= np.min(df.values)
    df = df[np.std(df, axis=1) / np.mean(np.std(df.values)) > 0.01]
        
    isEC = df.columns.get_level_values('celltype').values == 'EC'
    df = df.droplevel('celltype', axis=1).astype(float)
    df_EC = df[df.columns[isEC]]
    df_other = df[df.columns[~isEC]]

    nBatches = 10
    EC_batches = np.hstack([v + '_' + str(i) for i, v in enumerate(np.array_split(df_EC.columns.get_level_values('batch').values, nBatches))])
    other_batches = np.hstack([v + '_' + str(i) for i, v in enumerate(np.array_split(df_other.columns.get_level_values('batch').values, nBatches))])
    np.random.shuffle(EC_batches)
    np.random.shuffle(other_batches)

    df_EC.columns = pd.MultiIndex.from_arrays([EC_batches, df_EC.columns.get_level_values('cell')], names=['batch', 'cell'])
    df_other.columns = pd.MultiIndex.from_arrays([other_batches, df_other.columns.get_level_values('cell')], names=['batch', 'cell'])

    print(df_EC)
    print(df_other)

    return df_EC, df_other

# End: Part specific for each analyses ------------------------------------------------------------


# Begin: Part not to be edited for any specific analyses ------------------------------------------

from analysisPipeline import analyze, compareTwoCases

if not os.path.exists(workingDir):
    os.makedirs(workingDir)

byBatchesDir = workingDir + 'byBatches/'
bootstrapDir = workingDir + 'bootstrap/'
dataSaveName = workingDir + 'data.h5'
dendroDataName = workingDir + 'bootstrap_experiments_dendro_data.h5'

bootstrapExperiments = list(range(0, 100))

def prepareDEG(dfa, dfb, saveName):
    
    print('Saving expression data', flush=True)
    dfa.to_hdf(saveName, key='df', mode='a', complevel=4, complib='zlib')

    genes = []
    batches = []
    for batch in np.unique(dfa.columns.get_level_values('batch').values):
        df_ttest = pd.DataFrame(index=dfa.index, columns=['statistic', 'pvalue'])
        ttest = scipy.stats.ttest_ind(dfa.xs(batch, level='batch', axis=1, drop_level=False).values, 
                                        dfb.xs(batch, level='batch', axis=1, drop_level=False).values, axis=1)
        df_ttest['statistic'] = ttest[0]
        df_ttest['pvalue'] = ttest[1]
        df_ttest = df_ttest.sort_values('statistic', ascending=False).dropna()
        genes.append(df_ttest.loc[df_ttest['pvalue'] <= 10**-3]['statistic'].index.values)
        batches.append(batch)

    ugenes = []
    for i, v in enumerate(genes):
        ugenes.extend(v)
    ugenes = np.unique(ugenes)
    print('\nBatches:', len(genes), 'Unique genes:', len(ugenes))

    df = pd.DataFrame(index=range(len(ugenes)), columns=batches)
    for i, v in enumerate(genes):
        df.iloc[:len(v), i] = v

    df_ranks = pd.DataFrame(index=ugenes, columns=batches)
    for i, v in enumerate(genes):
        df_ranks.iloc[:, i] = pd.Index(df.iloc[:, i]).reindex(df_ranks.index)[1]

    print('Saving ranks data', flush=True)
    df_ranks.to_hdf(saveName, key='df_ranks', mode='a', complevel=4, complib='zlib')
    print(df_ranks)

    return

def prepareBootstrapExperiments(sourceDir, saveDir, ids = [], majorMetric = 'correlation', allDataToo = True, df_ranks = None):

    try:
        print('Reading precalculated data from %s' % sourceDir, flush=True)
        df_measure = pd.read_hdf(os.path.join(sourceDir, 'metricsFile.h5'), key=majorMetric)
        print('df_measure', '%1.1fMb'%(sys.getsizeof(df_measure) / 1024**2), flush=True)

        df_fraction =  pd.read_hdf(os.path.join(sourceDir, 'perGeneStats.h5'), key='df_fraction')
        df_median_expr =  pd.read_hdf(os.path.join(sourceDir, 'perGeneStats.h5'), key='df_expression')
        se_count =  pd.read_hdf(os.path.join(sourceDir, 'perGeneStats.h5'), key='se_count')

        saveSubDirs = ['Experiment %s' % (id + 1) for id in ids]
        if allDataToo:
            saveSubDirs = ['All'] + saveSubDirs

        for saveSubDir in saveSubDirs:
            try:
                print('\n', saveSubDir, flush=True)
                if not os.path.exists(os.path.join(saveDir, saveSubDir)):
                    os.makedirs(os.path.join(saveDir, saveSubDir))

                if saveSubDir == 'All':
                    batches = df_measure.columns
                else:
                    batches = np.random.choice(df_measure.columns, size=len(df_measure.columns), replace=True)

                np.savetxt(os.path.join(saveDir, saveSubDir, 'batches.txt'), batches, fmt='%s')

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
        
                np.savetxt(os.path.join(saveDir, saveSubDir, 'size.txt'), [df_fraction_temp.shape[0], se_count_temp.sum()], fmt='%i')

                if not df_ranks is None:
                    df_ranks_temp = df_ranks[batches]
                    df_ranks_temp.columns = df_ranks_temp.columns + '_' + np.array(range(len(df_ranks_temp.columns))).astype(str)
                    df_ranks_temp = df_ranks_temp.median(axis=1).sort_values()
                    df_ranks_temp.to_hdf(os.path.join(saveDir, saveSubDir, 'perGeneStats.h5'), key='df_ranks', mode='a', complevel=4, complib='zlib')

            except Exception as exception:
                print(exception)

    except Exception as exception:
        print(exception)

    return

def runPairOfExperiments(args):

    saveDir, saveSubDir = args

    print(saveDir, saveSubDir, flush=True)

    try:
        comparisonName = os.path.join(saveDir, saveSubDir, 'comparison')

        analyze(None, receptorsListHugo_2555, gEC23, [], 'correlation', toggleAdjustText=False, noPlot=True, panels=[],
                suffix=saveSubDir, saveDir=os.path.join(saveDir, saveSubDir), printStages=False,
                toggleCalculateMajorMetric=False, toggleExportFigureData=True, toggleCalculateMeasures=True)

        compareTwoCases(os.path.join(saveDir, saveSubDir, ''), otherCaseDir, name1='name1', name2='name2', saveName=comparisonName)

        additionalData = externalPanelsData.copy()
        additionalData.update({'conservedGenes': pd.read_excel(comparisonName + '.xlsx', index_col=1, header=0)['Inter-measures.T50_common_count']})

        analyze(None, receptorsListHugo_2555, gEC23, [], 'correlation', toggleAdjustText=False, dpi=300,
                suffix=saveSubDir, saveDir=os.path.join(saveDir, saveSubDir), printStages=False,
                toggleCalculateMajorMetric=False, toggleExportFigureData=True, toggleCalculateMeasures=False, externalPanelsData=additionalData) 

    except Exception as exception:
        print(exception)

    return

def analyzeCombinationVariant(variant):

    def getPeaksLists(df_temp):
        peaksListsMerged = []
        peaksLists = []
        experiments = []
        for experiment in np.unique(df_temp.index.get_level_values('experiment')):
            se = df_temp.xs(experiment, level='experiment', axis=0)
            genesInHalfPeak = getGenesOfPeak(se)

            peaksListsMerged.extend(genesInHalfPeak)
            peaksLists.append(genesInHalfPeak)
            experiments.append(experiment)

        return peaksListsMerged, peaksLists, experiments

    print('Variant:', variant)

    df = pd.read_hdf(dendroDataName, key='df').fillna(0).set_index('gene', append=True).droplevel('order')[variant]
            
    variabilitySavePath = extractBootstrapVariability(variant, filePath=dendroDataName,
                        savePath=workingDir + '/variability_' + variant.replace(' ', '-') + '.xlsx')

    writer = pd.ExcelWriter(workingDir + '%s bootstrap_in-peak_genes_SD.xlsx' % variant)

    result = dict()
    for species in ['species']:
        try:
            listsMerged, lists, experiments = getPeaksLists(df.xs(species, level='species', axis=0).copy())
            listsSizes = [len(item) for item in lists]
            #print('\n', species, 'mean, median:', np.mean(listsSizes), np.median(listsSizes))

            umL = np.unique(listsMerged, return_counts=True)
            se = pd.Series(index=umL[0], data=umL[1]/len(experiments)).sort_values(ascending=False)

            filePath = os.path.join(bootstrapDir + '/All/', 'dendrogram-heatmap-correlation-data.xlsx')
            peakGenesAll = getGenesOfPeak(pd.read_excel(filePath, index_col=0, header=0, sheet_name='Cluster index')[variant])
            df_res = pd.concat([se, se], axis=1, sort=False)
            df_res.iloc[:, 1] = np.where(np.isin(df_res.index.values, peakGenesAll), 1, np.nan)
            df_res.columns = ['Bootstrap', 'In all']

            try:
                dfv = pd.read_excel(variabilitySavePath, sheet_name=species, header=0, index_col=0).reindex(df_res.index)
                df_res = pd.concat([df_res, dfv], axis=1, sort=False)
            except Exception as exception:
                print('ERROR reading variability file:', exception)

            df_res.to_excel(writer, species)
            result.update({(species, variant): df_res})

            df_exp = pd.DataFrame(index=range(1000), columns=experiments)
            for i, col in enumerate(df_exp.columns):
                df_exp.iloc[:len(lists[i]), i] = lists[i]
            df_exp.to_excel(writer, 'Bootstrap lists ' + species, index=False)

        except Exception as exception:
            print('ERROR:', exception)

    writer.save()

    return result

# End: Part not to be edited for any specific analyses --------------------------------------------

if __name__ == '__main__':

    # Step 0. Preapre gene expression data and DEG ranking
    if True:
        df_EC, df_other = prepareInputData()

        prepareDEG(df_EC, df_other, dataSaveName)

    # Step 1. Generate per-batch correlations
    if True:
        analyze(pd.read_hdf(dataSaveName, key='df'), 
                receptorsListHugo_2555, gEC23, [], 'correlation', exprCutoff=0.02, suffix='all', saveDir=byBatchesDir,
                toggleCalculateMajorMetric=True, toggleExportFigureData=True, toggleCalculateMeasures=True, toggleGroupBatches=False, nCPUs=nCPUs)

    # Step 2. Generate data for bootstrap experiments
    if True:
        np.random.seed(0)
        prepareBootstrapExperiments(byBatchesDir, bootstrapDir, ids=bootstrapExperiments, allDataToo=True, df_ranks=pd.read_hdf(dataSaveName, key='df_ranks'))

    # Step 3. Analyze bootstrap experiments and collect all dendro data
    if True:
        saveSubDirs = ['All'] + ['Experiment %s' % (id + 1) for id in bootstrapExperiments]
        pool = multiprocessing.Pool(processes=nCPUs)
        pool.map(runPairOfExperiments, [(bootstrapDir, saveSubDir) for saveSubDir in saveSubDirs])
        pool.close()
        pool.join()

        dfs = []
        for id in bootstrapExperiments:
            saveSubDir = 'Experiment %s' % (id + 1)

            filePath = os.path.join(bootstrapDir, saveSubDir, 'dendrogram-heatmap-correlation-data.h5')
            try:
                df_temp = pd.read_hdf(filePath, key='df_C')
                df_temp.index.name = 'gene'
                df_temp = df_temp.reset_index()
                df_temp = pd.concat([df_temp], keys=[saveSubDir], axis=0, sort=False)
                df_temp = pd.concat([df_temp], keys=['species'], axis=0, sort=False)
                df_temp.index.names = ['species', 'experiment', 'order']
                dfs.append(df_temp)

            except Exception as exception:
                print(exception)
                pass

        dfs = pd.concat(dfs, axis=0, sort=False)
        print(dfs)
        dfs.to_hdf(dendroDataName, key='df', mode='a', complevel=4, complib='zlib')

    # Step 4. Anlyze peaks from bootstrap experiments
    if True:
        totalResults = dict()
        for variant in ['Avg combo4avgs', 'Avg combo3avgs']: 
            totalResults.update(analyzeCombinationVariant(variant))

        df_totalResults = []
        for species, variant in totalResults.keys():
            print(variant, species)
            df_temp = totalResults[(species, variant)]
            df_temp.index.name = 'genes'
            df_temp = df_temp.reset_index()
            df_temp = pd.concat([df_temp], keys=[variant], names=['variant'], axis=1)
            df_temp = pd.concat([df_temp], keys=[species], names=['species'], axis=1)
            df_totalResults.append(df_temp)

        df_totalResults = pd.concat(df_totalResults, axis=1, sort=False)
        df_totalResults.to_excel(workingDir + 'variants.xlsx')
        print(df_totalResults)

        df_totalResults = df_totalResults[df_totalResults.columns[np.isin(df_totalResults.columns.get_level_values(-1).values, ['genes', 'Bootstrap'])]]
        df_totalResults.to_excel(workingDir + 'variants_cut.xlsx')
        print(df_totalResults)

    # Step 5. Dendrogram randomization
    if True:
        df = pd.read_excel(os.path.join(bootstrapDir + 'All/dendrogram-heatmap-correlation-data.xlsx'), index_col=0, header=0, sheet_name='Cluster index')
            
        scramble(df.copy(), ['Markers', 'Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], workingDir=workingDir + 'random/combo4/', M=20)
        scramble(df.copy(), ['Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], workingDir=workingDir + 'random/combo3/', M=20)

    # Step 6. Make "All" plot with text adjusted
    if True:
        additionalData = externalPanelsData.copy()
        additionalData.update({'conservedGenes': pd.read_excel(os.path.join(bootstrapDir, 'All', 'comparison.xlsx'), index_col=1, header=0)['Inter-measures.T50_common_count']})

        analyze(None, receptorsListHugo_2555, gEC23, [], 'correlation', toggleAdjustText=True, dpi=300,
                suffix='All', saveDir=os.path.join(bootstrapDir, 'All/'), printStages=True,
                toggleCalculateMajorMetric=False, toggleExportFigureData=True, toggleCalculateMeasures=False, externalPanelsData=additionalData) 