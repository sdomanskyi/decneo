from commonFunctions import *
from analysisPipeline import analyze, compareTwoCases


workingDir = 'results/Endothelial by PanglaoDB definition/'

def collect_EC():

    dataFile = workingDir + 'PanglaoDB_EC.h5'
    dataPreFile = workingDir + 'PanglaoDB_EC_temp.h5'

    EC_cells_from_Alex = pd.read_csv(workingDir + 'Alex_Endothelial_Cells.csv', index_col=0, header=0, dtype=str).stack().dropna().reset_index().set_index(['level_1', 0]) # 87 SRA, 584 SRS, 138141 cells
    EC_cells_from_Alex.index.names = ['batch', 'cell']
    
    batches = np.unique(EC_cells_from_Alex.index.get_level_values('batch').values)
    for ibatch, batch in enumerate(batches):
        if not KeyInStore(batch, dataPreFile):
            try:
                print('Processing:', batch, '\t', ibatch, 'of', len(batches), flush=True)

                df_expr = pd.read_hdf(os.path.join(RDataDirName, batch + '.sparse.RData.h5'), key='df').reindex(EC_cells_from_Alex.xs(batch, level='batch').index, axis=1)
                df_expr.index = pd.Series(df_expr.index.values).replace(Mouse_to_Human_HUGO_conversion).values

                df_expr = df_expr.loc[~df_expr.index.duplicated(keep='first')]
                df_expr = df_expr.T.loc[~df_expr.T.index.duplicated(keep='first')].T
                df_expr = df_expr[df_expr.sum(axis=1) > 0].fillna(0).astype(int)

                df_expr = pd.concat([df_expr], keys=[batch], names=['batch'], axis=1, sort=True)
                df_expr.columns.names = ['batch', 'cell']

                print(df_expr, flush=True)

                df_expr.to_hdf(dataPreFile, key=batch, mode='a', complevel=4, complib='zlib')

            except Exception as exception:
                print('Error while processing file:\n', exception)

    keys = KeysOfStore(dataPreFile)
    print('Reading %s keys' % (len(keys)), flush=True)

    df_all = []
    for key in keys:
        print('\tReading:', key, flush=True)
        df_all.append(pd.read_hdf(dataPreFile, key=key))

    print('Merging data', flush=True)
    df_all = pd.concat(df_all, axis=1).fillna(0.).astype(int)

    print('Recording all samples data to hdf', flush=True)
    df_all.to_hdf(dataFile, key='df', mode='a', complevel=4, complib='zlib')

    return

def normalize_EC():

    for species in ['Homo sapiens', 'Mus musculus'][1:]:
        batches = pd.read_excel(workingDir + 'batches.xlsx', index_col=0, header=0)['species']
        batches = batches.loc[batches==species].index.values
        print('All %s batches:' % species, len(batches))

        try:
            print('Loading data', flush=True)
            df = pd.read_hdf(workingDir + 'PanglaoDB_EC.h5')
            df = df[df.columns[np.isin(df.columns.get_level_values('batch').values, batches)]]
            print(df.shape)

            print('Converting index genes', flush=True)
            df.index = DigitalCellSorter.DigitalCellSorter().gnc.Convert(df.index.values.tolist(), 'alias', 'hugo', returnUnknownString=False)
            df = df.loc[~df.index.duplicated(keep='first')].astype(float)

            print('Scaling and log-transforming cells', flush=True)
            df /= df.sum(axis=0) * 0.0001
            df = np.log2(df.replace(0., np.min(df.values[df.values > 0.])))
            df -= np.min(df.values)

            print('Removing constant genes', flush=True)
            s0 = set(df.index.values)
            df = df[np.std(df, axis=1) / np.mean(np.std(df.values)) > 0.01]
            print(s0.difference(set(df.index.values)))
            print(df)

            df.to_hdf(workingDir + 'PanglaoDB_EC_normed_and_filtered.h5', key=species, mode='a', complevel=4, complib='zlib')
        except Exception as exception:
            print('ERROR:', exception, flush=True)
            print('\n\n', flush=True)

    return df

def ttest_EC():
    
    dataFile = workingDir + 'PanglaoDB_EC_ttest.h5'
    dataPreFile = workingDir + 'PanglaoDB_EC_ttest_temp.h5'

    if False:
        EC_cells_from_Alex = pd.read_csv(workingDir + 'Alex_Endothelial_Cells.csv', index_col=0, header=0, dtype=str).stack().dropna().reset_index().set_index(['level_1', 0]) # 87 SRA, 584 SRS, 138141 cells
        EC_cells_from_Alex.index.names = ['batch', 'cell']
    
        batches = np.unique(EC_cells_from_Alex.index.get_level_values('batch').values)
        for ibatch, batch in enumerate(batches):
            if not KeyInStore(batch, dataPreFile):
                if ibatch >= 0:
                    try:
                        print('Processing:', batch, '\t', ibatch, 'of', len(batches), flush=True)

                        df_expr = pd.read_hdf(os.path.join(RDataDirName, batch + '.sparse.RData.h5'), key='df')
                        df_expr.index = pd.Series(df_expr.index.values).replace(Mouse_to_Human_HUGO_conversion).values

                        cellEC = EC_cells_from_Alex.xs(batch, level='batch').index
                        cells = pd.Index(np.random.choice(df_expr.columns.difference(cellEC), 1000, replace=False)).append(cellEC)
                        df_expr = df_expr.reindex(cells, axis=1)

                        df_expr = df_expr.loc[~df_expr.index.duplicated(keep='first')]
                        df_expr = df_expr.T.loc[~df_expr.T.index.duplicated(keep='first')].T
                        df_expr = df_expr[df_expr.sum(axis=1) > 0].fillna(0).astype(int)
                        df_expr = pd.concat([df_expr], keys=[batch], names=['batch'], axis=1, sort=True)
                        df_expr.columns.names = ['batch', 'cell']
                
                        df_EC = df_expr.reindex(EC_cells_from_Alex.xs(batch, level='batch', drop_level=False).index, axis=1)
                        df_other = df_expr[df_expr.columns.difference(df_EC.columns)]

                        print(df_expr.shape, df_EC.shape, df_other.shape, flush=True)

                        df_ttest = pd.DataFrame(index=df_EC.index, columns=['statistic', 'pvalue'])
                        ttest = scipy.stats.ttest_ind(df_EC.values, df_other.values, axis=1)
                        df_ttest['statistic'] = ttest[0]
                        df_ttest['pvalue'] = ttest[1]
                        df_ttest = df_ttest.sort_values('statistic', ascending=False).dropna()

                        df_ttest.to_hdf(dataPreFile, key=batch, mode='a', complevel=4, complib='zlib')

                    except Exception as exception:
                        print('Error while processing file:\n', exception)

    if False:
        keys = KeysOfStore(dataPreFile)
        print('Reading %s keys' % (len(keys)), flush=True)

        df_all = []
        for batch in keys:
            print('\tReading:', batch, flush=True)
            df_temp = pd.read_hdf(dataPreFile, key=batch)
            df_temp = pd.concat([df_temp], keys=[batch[1:]], names=['batch'], axis=1, sort=True)
            df_all.append(df_temp)

        print('Merging data', flush=True)
        df_all = pd.concat(df_all, axis=1).fillna(0.).astype(int)

        print('Recording all samples data to hdf', flush=True)
        df_all.to_hdf(dataFile, key='df', mode='a', complevel=4, complib='zlib')

    if False:
        batches = pd.read_excel(workingDir + 'batches.xlsx', index_col=0, header=0)['species']
        humanBatches = batches.loc[batches=='Homo sapiens'].index.values
        mouseBatches = batches.loc[batches=='Mus musculus'].index.values
        
        for species in ['Homo sapiens', 'Mus musculus']:
            batches = humanBatches if species == 'Homo sapiens' else mouseBatches
            print(species, len(batches))

            fileName = workingDir + 'PanglaoDB_EC_ttest_temp.h5'

            genes = []
            genes_batches = []
            for i, batch in enumerate(batches):
                try:
                    df_ttest = pd.read_hdf(fileName, key=batch)
                    genes.append(df_ttest.loc[df_ttest['pvalue'] <= 10**-3]['statistic'].index.values)
                    genes_batches.append(batch)
                except Exception as exception:
                    pass

            ugenes = []
            for i, v in enumerate(genes):
                ugenes.extend(v)
            ugenes = np.unique(ugenes)
            print(len(genes), len(ugenes))

            df = pd.DataFrame(index=range(len(ugenes)), columns=genes_batches)
            for i, v in enumerate(genes):
                df.iloc[:len(v), i] = v

            df_ranks = pd.DataFrame(index=ugenes, columns=genes_batches)
            for i, v in enumerate(genes):
                df_ranks.iloc[:, i] = pd.Index(df.iloc[:, i]).reindex(df_ranks.index)[1]

            print(df_ranks)
            df_ranks.to_hdf(workingDir + 'PanglaoDB_ttest_ranks_per_batch %s.h5' % species, key='df', mode='a', complevel=4, complib='zlib')

            df_ranks = df_ranks.apply(np.median, axis=1).sort_values()
            df_ranks = df_ranks[df_ranks > -1][:]
            print(df_ranks)

            df_ranks.to_excel(workingDir + 'from_ranks_Panglao %s.xlsx' % species)

    return

def generatePerBatchCorrelation_analyze_EC():

    fname = workingDir + 'PanglaoDB_EC_normed_and_filtered.h5'

    for species in ['Homo sapiens', 'Mus musculus']:
        print('Analyzing all cells of %s:' % species, flush=True)

        try:
            analyze(pd.read_hdf(fname, key=species), receptorsListHugo_2555, gECs, gECi, 'correlation',
                    suffix='All, ' + species, saveDir=workingDir + 'byBatches/%s/' % species,
                    toggleCalculateMajorMetric=True, toggleExportFigureData=True, toggleCalculateMeasures=True, toggleGroupBatches=False,
                    diffExpressedGenes=None, conservation=None, toggleAdjustText=False, noPlot=True, panels=[])

        except Exception as exception:
            print('ERROR:', exception, flush=True)
            print('\n\n', flush=True)

    return

def prepareBootstrapExperiments2(sourceDir, saveDir, ids = [], majorMetric = 'correlation', allDataToo = False, df_ranks = None):

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
        
                if not df_ranks is None:
                    df_ranks_temp = df_ranks[batches]
                    df_ranks_temp.columns = df_ranks_temp.columns + '_' + np.array(range(len(df_ranks_temp.columns))).astype(str)
                    df_ranks_temp = df_ranks_temp.median(axis=1).sort_values()
                    df_ranks_temp.to_hdf(os.path.join(saveDir, saveSubDir, 'perGeneStats.h5'), key='df_ranks', mode='a', complevel=4, complib='zlib')

                np.savetxt(os.path.join(saveDir, saveSubDir, 'size.txt'), [df_fraction_temp.shape[0], se_count_temp.sum()], fmt='%i')

            except Exception as exception:
                print(exception)

    except Exception as exception:
        print(exception)

    return

def prepareBootstrapExperiments3(sourceDir, saveDir, ids = [], majorMetric = 'correlation', allDataToo = False, df_ranks = None):

    saveSubDirs = ['Experiment %s' % (id + 1) for id in ids]
    if allDataToo:
        saveSubDirs = ['All'] + saveSubDirs

    for saveSubDir in saveSubDirs:
        try:
            print('\n', saveSubDir, flush=True)
            batches = np.loadtxt(os.path.join(saveDir, saveSubDir, 'batches.txt'), dtype=str)
            df_ranks_temp = df_ranks[batches]
            df_ranks_temp.columns = df_ranks_temp.columns + '_' + np.array(range(len(df_ranks_temp.columns))).astype(str)
            df_ranks_temp = df_ranks_temp.median(axis=1).sort_values()
            df_ranks_temp.to_hdf(os.path.join(saveDir, saveSubDir, 'perGeneStats.h5'), key='df_ranks', mode='a', complevel=4, complib='zlib')

        except Exception as exception:
            print(exception)

    return

def runPairOfExperiments(args):

    saveDir, saveSubDir = args
    print(saveDir, saveSubDir, flush=True)

    try:
        comparisonName = os.path.join(saveDir, 'Mus musculus', saveSubDir, 'comparison')

        if True:
            for species in ['Homo sapiens', 'Mus musculus']:
                #print('Processing:', species)
                analyze(None, receptorsListHugo_2555, gECs, gECi, 'correlation', toggleAdjustText=False, noPlot=True, panels=[],
                        suffix=saveSubDir + ', ' + species, saveDir=os.path.join(saveDir, species, saveSubDir), printStages=False,
                        toggleCalculateMajorMetric=False, toggleExportFigureData=True, toggleCalculateMeasures=True,
                        diffExpressedGenes=None, conservation=None)

        if True:
            compareTwoCases(os.path.join(saveDir, 'Homo sapiens', saveSubDir, ''), 
                            os.path.join(saveDir, 'Mus musculus', saveSubDir, ''), 
                            name1='human', name2='mouse', saveName=comparisonName)

            conservation = {'conservedGenes': pd.read_excel(comparisonName + '.xlsx', index_col=1, header=0)['Inter-measures.T50_common_count'],
                           'conservedMarkers': pd.read_excel(comparisonName + '.xlsx', index_col=1, header=0)['Inter-measures.EC23T50_common_count']}
    
            for species in ['Homo sapiens', 'Mus musculus']:
                #print('Re-processing:', species)
                analyze(None, receptorsListHugo_2555, gECs, gECi, 'correlation', toggleAdjustText=False, dpi=300,
                        suffix=saveSubDir + ', ' + species, saveDir=os.path.join(saveDir, species, saveSubDir), printStages=False,
                        toggleCalculateMajorMetric=False, toggleExportFigureData=True, toggleCalculateMeasures=False,
                        diffExpressedGenes=None, conservation=conservation) 

    except Exception as exception:
        print(exception)

    return


if __name__ == '__main__':

    if False:
        df = pd.read_excel('dev/PanglaoDBdata/df_cell_type_annotations.xlsx', index_col=[0, 1], header=0)[['Cell type annotation', 'Species']].reset_index()[['Species', 'Cell type annotation', 'SRA accession']].drop_duplicates()
        df = df.set_index(['Species', 'Cell type annotation'])['SRA accession']
        print(df)

        dfg = df.groupby(level=['Species', 'Cell type annotation'], axis=0).count().unstack('Species').sort_values('Mus musculus', ascending=False).reset_index()
        print(dfg)
        dfg.to_excel('dfgSRA.xlsx')
        exit()

    if False:
        df = pd.read_excel('Alona_DCS_all_clusters.xlsx', index_col=None, header=0)

        def case(df, i, j):

            se = df[df.columns[[i, j]]]
            se = se.set_index(se.columns[0])[se.columns[1]].dropna().astype(int)
            print(se)

            return se

        se_mA = case(df, 0, 1)
        se_hA = case(df, 2, 3)
        se_mD = case(df, 4, 5)
        se_hD = case(df, 6, 7)

        def pair(se1, se2):

            df = pd.DataFrame(index=pd.MultiIndex.from_arrays(np.unique(se1.values, return_counts=True)), 
                              columns=pd.MultiIndex.from_arrays(np.unique(se2.values, return_counts=True)))

            for i in df.index:
                for j in df.columns:
                    df.loc[i, j] = len(se1[se1==i[0]].index.intersection(se2[se2==j[0]].index))

            print(df)

            return df

        writer = pd.ExcelWriter('cross2.xlsx')
        pair(se_mA, se_hA).to_excel(writer, 'mA_hA', merge_cells=False)
        pair(se_mD, se_hD).to_excel(writer, 'mD_hD', merge_cells=False)
        pair(se_mA, se_mD).to_excel(writer, 'mA_mD', merge_cells=False)
        pair(se_hA, se_hD).to_excel(writer, 'hA_hD', merge_cells=False)
        writer.save()

        exit()

    nCPUs = 4 if platform.system() == "Windows" else 15
    bootstrapExperiments = list(range(0, 100))
    bootstrapDir = workingDir + 'EC all bootstrap 100 w21/'

    # Step 1. Generate data for bootstrap experiments
    if False:
        np.random.seed(0)

        try:
            batchID = eval(sys.argv[1]) # 0, 1, 2, 3
            np.random.seed(11 * batchID)
            bootstrapExperiments = bootstrapExperiments[12+22*batchID:12+22*(batchID+1)]
            print(bootstrapExperiments, flush=True)
        except:
            pass
        
        for species in ['Homo sapiens', 'Mus musculus'][1:]:
            df_ranks = pd.read_hdf(workingDir + 'PanglaoDB_ttest_ranks_per_batch %s.h5' % species, key='df')
            prepareBootstrapExperiments2(workingDir + 'byBatches/%s/' % species, bootstrapDir + '%s/' % species, 
                                         allDataToo=False, ids=bootstrapExperiments, df_ranks=df_ranks)

    # Step 2. Analyze bootstrap experiments
    if False:
        saveSubDirs = ['All'] + ['Experiment %s' % (id + 1) for id in bootstrapExperiments]
        pool = multiprocessing.Pool(processes=nCPUs)
        pool.map(runPairOfExperiments, [(bootstrapDir, saveSubDir) for saveSubDir in saveSubDirs])
        pool.close()
        pool.join()

    # Step 3. Collect all dendro data from 100 experiments
    if False:
        dfs = []
        for id in bootstrapExperiments:
            saveSubDir = 'Experiment %s' % (id + 1)
            print(saveSubDir)
            for species in ['Homo sapiens', 'Mus musculus']:
                filePath = os.path.join(bootstrapDir, species, saveSubDir, 'dendrogram-heatmap-correlation-data.xlsx')
                try:
                    df_temp = pd.read_excel(filePath, index_col=0, header=0, sheet_name='Cluster index')
                    df_temp.index.name = 'gene'
                    df_temp = df_temp.reset_index()
                    df_temp = pd.concat([df_temp], keys=[saveSubDir], axis=0, sort=False)
                    df_temp = pd.concat([df_temp], keys=[species], axis=0, sort=False)
                    df_temp.index.names = ['species', 'experiment', 'order']
                    dfs.append(df_temp)
                    print(df_temp)

                except Exception as exception:
                    print(exception)
                    pass

        dfs = pd.concat(dfs, axis=0, sort=False)
        print(dfs)
        dfs.to_hdf(bootstrapDir + 'bootstrap_100experiments_dendro_data.h5', key='df', mode='a', complevel=4, complib='zlib')

    # Step 4.A. Anlyze peaks dendro data from 100 experiments
    if True:
        def getGenesOfPeak(se):

            se /= se.max()
            i0 = np.argmax(se.values)
            genesInHalfPeak = [se.index[i0]]

            for i in range(i0 + 1, i0 + 10**3, 1):
                try:
                    if se.iloc[i] >= 0.5:
                        genesInHalfPeak.append(se.index[i])
                    else:
                        break
                except:
                    break

            for i in range(i0-1, 0, -1):
                try:
                    if se.iloc[i] >= 0.5:
                        genesInHalfPeak.append(se.index[i])
                    else:
                        break
                except:
                    break

            return genesInHalfPeak

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

        def runCol(variant):
            print('Variant:', variant)

            df = pd.read_hdf(bootstrapDir + 'bootstrap_100experiments_dendro_data.h5', key='df').fillna(0).set_index('gene', append=True).droplevel('order')[variant]

            writer = pd.ExcelWriter(bootstrapDir + '%s 100_bootstrap_in-peak_genes_SD.xlsx' % variant)

            for species in ['Mus musculus', 'Homo sapiens']:

                print('\n', species)

                listsMerged, lists, experiments = getPeaksLists(df.xs(species, level='species', axis=0).copy())
                listsSizes = [len(item) for item in lists]
                print('mean, median:', np.mean(listsSizes), np.median(listsSizes))

                umL = np.unique(listsMerged, return_counts=True)
                se = pd.Series(index=umL[0], data=umL[1]).sort_values(ascending=False)

                filePath = os.path.join(bootstrapDir + species + '/All/', 'dendrogram-heatmap-correlation-data.xlsx')
                peakGenesAll = getGenesOfPeak(pd.read_excel(filePath, index_col=0, header=0, sheet_name='Cluster index')[variant])
                se = pd.concat([se, se], axis=1, sort=False)
                se.iloc[:, 1] = np.where(np.isin(se.index.values, peakGenesAll), 1, np.nan)
                se.columns = ['Bootstrap', 'In all']
                se.to_excel(writer, species)

                df_exp = pd.DataFrame(index=range(1000), columns=experiments)
                for i, col in enumerate(df_exp.columns):
                    df_exp.iloc[:len(lists[i]), i] = lists[i]
                df_exp.to_excel(writer, 'Bootstrap lists ' + species, index=False)

            writer.save()

            return

        #runCol('Avg 13 Combination 3 avgs')
        runCol('Avg 11 Combination 4 avgs')
