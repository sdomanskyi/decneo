from commonFunctions import *
from analysisPipeline import analyze, compareTwoCases

if platform.system() == "Windows":
    workingDir = 'choroid Voigt bootstrap Enriched original/'
    bootstrapDir = workingDir + 'Panglao mouse vs Human choroid Voigt vs all/'
    fname_PanglaoDB_EC_byDCS = 'DCS output/PanglaoDB_EC_byDCS_normed_and_filtered_5k.h5'
    fname_Choroid_Voigt_EC_byDCS = 'd:/Projects/A_Endothelial/Voigt_choroid_EC_v2.h5'
    fname_Choroid_Voigt_AllCells_byDCS = 'd:/Projects/A_Endothelial/Eye/Voigt_Choroid.h5'
    fname_Voigt_Choroid_original_data_AH_subset = 'd:/Projects/A_Endothelial/Eye/Voigt_Choroid_original_data_AH_subset.h5'
    fname_Voigt_Choroid_original_data_AH_subset_ranks = 'd:/Projects/A_Endothelial/Eye/Voigt_Choroid_original_data_AH_subset_df_ranks.h5'
else:
    workingDir = 'choroid Voigt bootstrap Enriched original/'
    bootstrapDir = workingDir + 'Panglao mouse vs Human choroid Voigt vs all/'
    fname_PanglaoDB_EC_byDCS = 'DCS output/PanglaoDB_EC_byDCS_normed_and_filtered.h5'
    fname_Choroid_Voigt_EC_byDCS = '/mnt/research/piermarolab/data/Eye/Voigt_choroid_EC_v2.h5'
    fname_Choroid_Voigt_AllCells_byDCS = '/mnt/research/piermarolab/data/Eye/Voigt_Choroid.h5'
    fname_Voigt_Choroid_original_data_AH_subset = '/mnt/research/piermarolab/data/Eye/Voigt_Choroid_original_data_AH_subset.h5'
    fname_Voigt_Choroid_original_data_AH_subset_ranks = '/mnt/research/piermarolab/data/Eye/Voigt_Choroid_original_data_AH_subset_df_ranks.h5'

def preapreExpressionAndDEG():
    
    se_AH_labels = pd.read_csv('choroid compare labels to AH labels/Choroid_Cell_Types_by_AH.csv', index_col=0, header=0)['x']
    se_AH_labels = se_AH_labels[se_AH_labels == 'Endothelial']
    se_AH_labels[:] = pd.Series(data=se_AH_labels.index.str.split('_', expand=True).droplevel(3).droplevel(3).values).apply('_'.join)
    se_AH_labels.index = pd.Series(data=se_AH_labels.index.str.split('_', expand=True).droplevel(0).droplevel(0).droplevel(0).values).apply('_'.join)
    se_AH_labels.index = pd.MultiIndex.from_arrays([se_AH_labels.values, se_AH_labels.index], names=['batch', 'cell'])
    se_AH_labels[:] = 'Endothelial'
    se_AH_labels.name = 'celltype'
    se_AH_labels = se_AH_labels.to_frame().set_index(['celltype'], append=True)
    cells = se_AH_labels.index.droplevel('celltype')
    print(len(cells), cells[0], cells.names)

    df_all = pd.read_hdf('d:/Projects/A_Endothelial/Eye/GSE135922 Voigt 2020 choroid/VoigtChoroid_14samples_original.h5', key='df')
    df_all.columns = df_all.columns.droplevel('cluster').droplevel('celltype')

    df_EC = df_all[cells].fillna(0.)
    df_EC = df_EC.loc[df_EC.sum(axis=1) > 0.]
    print(df_EC)

    df_other = df_all[df_all.columns.difference(df_EC.columns)].reindex(df_EC.index)
    print(df_other)

    genes = []
    batches = []
    for batch in np.unique(df_EC.columns.get_level_values('batch').values):
        df_ttest = pd.DataFrame(index=df_EC.index, columns=['statistic', 'pvalue'])
        ttest = scipy.stats.ttest_ind(df_EC.xs(batch, level='batch', axis=1, drop_level=False).values, 
                                        df_other.xs(batch, level='batch', axis=1, drop_level=False).values, axis=1)
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
    df_ranks.to_hdf(fname_Voigt_Choroid_original_data_AH_subset_ranks, key='df', mode='a', complevel=4, complib='zlib')
    print(df_ranks)

    print('Saving expression data', flush=True)
    df.to_hdf(fname_Voigt_Choroid_original_data_AH_subset, key='df', mode='a', complevel=4, complib='zlib')
    print(df)

    return

def prepareBootstrapExperiments(sourceDir, saveDir, ids = [], majorMetric = 'correlation', allDataToo = False, df_ranks = None):

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
        comparisonName = os.path.join(saveDir, 'Homo sapiens', saveSubDir, 'comparison')

        if True:
            species = 'Homo sapiens'
            analyze(None, receptorsListHugo_2555, gEC23, [], 'correlation', toggleAdjustText=False, noPlot=True, panels=[],
                    suffix=saveSubDir + ', ' + species, saveDir=os.path.join(saveDir, species, saveSubDir), printStages=False,
                    toggleCalculateMajorMetric=False, toggleExportFigureData=True, toggleCalculateMeasures=True)

        if True:
            compareTwoCases(os.path.join(saveDir, 'Homo sapiens', saveSubDir, ''), 
                            'workflow/PanglaoDB EC all cells w21/Mus musculus/', 
                            name1='human', name2='mouse', saveName=comparisonName)

            additionalData = externalPanelsData.copy()
            additionalData.update({'conservedGenes': pd.read_excel(comparisonName + '.xlsx', index_col=1, header=0)['Inter-measures.T50_common_count'],
                                  'conservedMarkers': pd.read_excel(comparisonName + '.xlsx', index_col=1, header=0)['Inter-measures.EC23T50_common_count']})

            species = 'Homo sapiens'
            analyze(None, receptorsListHugo_2555, gEC23, [], 'correlation', toggleAdjustText=False, dpi=300,
                    suffix=saveSubDir + ', ' + species, saveDir=os.path.join(saveDir, species, saveSubDir), printStages=False,
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

    df = pd.read_hdf(bootstrapDir + 'bootstrap_100experiments_dendro_data.h5', key='df').fillna(0).set_index('gene', append=True).droplevel('order')[variant]
            
    variabilitySavePath = extractBootstrapVariability(variant, filePath=bootstrapDir + 'bootstrap_100experiments_dendro_data.h5',
                        savePath=bootstrapDir + '/variability_' + variant.replace(' ', '-') + '.xlsx')

    writer = pd.ExcelWriter(bootstrapDir + '%s 100_bootstrap_in-peak_genes_SD.xlsx' % variant)

    result = dict()
    for species in np.unique(df.index.get_level_values('species')):
        try:
            listsMerged, lists, experiments = getPeaksLists(df.xs(species, level='species', axis=0).copy())
            listsSizes = [len(item) for item in lists]
            #print('\n', species, 'mean, median:', np.mean(listsSizes), np.median(listsSizes))

            umL = np.unique(listsMerged, return_counts=True)
            se = pd.Series(index=umL[0], data=umL[1]/len(experiments)).sort_values(ascending=False)

            filePath = os.path.join(bootstrapDir + species + '/All/', 'dendrogram-heatmap-correlation-data.xlsx')
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

def analyzeVoigtChoroid():

    comparisonName = bootstrapDir + 'Homo sapiens/All/comparison'
    compareTwoCases(bootstrapDir + 'Homo sapiens/All/', 'workflow/PanglaoDB EC all cells w21/Mus musculus/', saveName=comparisonName)

    additionalData = externalPanelsData.copy()
    additionalData.update({'conservedGenes': pd.read_excel(comparisonName + '.xlsx', index_col=1, header=0)['Inter-measures.T50_common_count'],
                          'conservedMarkers': pd.read_excel(comparisonName + '.xlsx', index_col=1, header=0)['Inter-measures.EC23T50_common_count']})

    print('Re-plotting choroid Voigt')
    analyze(None, receptorsListHugo_2555, gECs, gECi, 'correlation', exprCutoff=0.01,
            suffix='Choroid_Homo sapiens', saveDir=bootstrapDir + 'Homo sapiens/All/', toggleIncludeHeatmap=True,
            toggleCalculateMajorMetric=False, toggleExportFigureData=True, toggleCalculateMeasures=False, externalPanelsData=additionalData)

    return


if __name__ == '__main__':

    nCPUs = 4 if platform.system() == "Windows" else 10
    bootstrapExperiments = list(range(0, 100))

    # Step 0. Preapre gene expression data and DEG ranking
    if False:

        preapreExpressionAndDEG()

    # Step 1. Generate per-batch correlations
    if False:
        analyze(pd.read_hdf(fname_Voigt_Choroid_original_data_AH_subset, key='df'), 
                receptorsListHugo_2555, gEC23, [], 'correlation', exprCutoff=0.01,
                suffix='Choroid_Voigt', saveDir=workingDir + 'ChoroidbyBatches/',
                toggleCalculateMajorMetric=True, toggleExportFigureData=True, toggleCalculateMeasures=True, toggleGroupBatches=False,
                toggleAdjustText=False, noPlot=True, panels=[])

    # Step 2. Generate data for bootstrap experiments
    if False:
        np.random.seed(5642)
        prepareBootstrapExperiments(workingDir + 'ChoroidbyBatches/', bootstrapDir + 'Homo sapiens/', ids=bootstrapExperiments, allDataToo=True, df_ranks=pd.read_hdf(fname_Voigt_Choroid_original_data_AH_subset_ranks, key='df'))

    # Step 3. Analyze bootstrap experiments and collect all dendro data
    if False:
        saveSubDirs = ['All'] + ['Experiment %s' % (id + 1) for id in bootstrapExperiments]
        pool = multiprocessing.Pool(processes=nCPUs)
        pool.map(runPairOfExperiments, [(bootstrapDir, saveSubDir) for saveSubDir in saveSubDirs])
        pool.close()
        pool.join()

        dfs = []
        for id in bootstrapExperiments[:]:
            saveSubDir = 'Experiment %s' % (id + 1)
            for species in ['Homo sapiens']:
                filePath = os.path.join(bootstrapDir, species, saveSubDir, 'dendrogram-heatmap-correlation-data.h5')
                try:
                    df_temp = pd.read_hdf(filePath, key='df_C')
                    df_temp.index.name = 'gene'
                    df_temp = df_temp.reset_index()
                    df_temp = pd.concat([df_temp], keys=[saveSubDir], axis=0, sort=False)
                    df_temp = pd.concat([df_temp], keys=[species], axis=0, sort=False)
                    df_temp.index.names = ['species', 'experiment', 'order']
                    dfs.append(df_temp)

                except Exception as exception:
                    print(exception)
                    pass

        dfs = pd.concat(dfs, axis=0, sort=False)
        print(dfs)
        dfs.to_hdf(bootstrapDir + 'bootstrap_100experiments_dendro_data.h5', key='df', mode='a', complevel=4, complib='zlib')

    # Step 4. Anlyze peaks from bootstrap experiments
    if False:
        totalResults = dict()
        #for variant in ['Avg combo7avgs', 'Avg combo6-3avgs', 'Avg combo6-2avgs', 'Avg combo6-1avgs', 'Avg combo5-3avgs', 'Avg combo5-2avgs', 'Avg combo5-1avgs', 'Avg combo4avgs', 'Avg combo3avgs', 'Avg combo2avgs','Avg Markers','Avg Binomial -log(pvalue)','Avg Top50 overlap','Avg Fraction'][:]:
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
        df_totalResults.to_excel(bootstrapDir + 'variants.xlsx')
        print(df_totalResults)

        df_totalResults = df_totalResults[df_totalResults.columns[np.isin(df_totalResults.columns.get_level_values(-1).values, ['genes', 'Bootstrap'])]]
        df_totalResults.to_excel(bootstrapDir + 'variants_cut.xlsx')
        print(df_totalResults)

    # Step 5A. Dendrogram randomization
    if False:
        df = pd.read_excel(os.path.join(bootstrapDir + 'Homo sapiens/All/dendrogram-heatmap-correlation-data.xlsx'), index_col=0, header=0, sheet_name='Cluster index')
            
        scramble(df.copy(), ['Markers', 'Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], workingDir=bootstrapDir + 'random/Homo sapiens/combo4/', M=20)
        scramble(df.copy(), ['Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], workingDir=bootstrapDir + 'random/Homo sapiens/combo3/', M=20)

    # For testing purposes
    if False:
        analyzeVoigtChoroid()
        runPairOfExperiments((bootstrapDir, 'All')); exit()
