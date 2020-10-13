''' File holding Analysis class 
'''

from .commonFunctions import *

class Analysis():

    '''Class of analysis and visualization functions for scRegulation
    
    Parameters:
        workingDir: str, Default ''
            Directory to retrieve and save file and results to
        
        otherCaseDir: str, Default ''
            Directory holding comparison data
        
        genesOfInterest: list, Default None
            Particular genes to analyze

        knownRegulators: list, Default None
            TEXT

        nCPUs: int, Default 1
            Number of CPUs to use for multiprocessing

        panels: list, Default None
            Particular measurements to analyze

        nBootstrap: int, Default 100
            Number of bootstraps to perform

        majorMetric: str, Default 'correlation'
            Metric name (e.g. ‘correlation’)

        perEachOtherCase: boolean, False 
            TEXT

        metricsFile: str, 'metricsFile.h5' 
            Name of file holding measurement data for specified metric

        seed: int, None 
            Used to set randomness deterministic
       
    '''

    def __init__(self, workingDir = '', otherCaseDir = '', genesOfInterest = None, knownRegulators = None, nCPUs = 1, panels = None, nBootstrap = 100, majorMetric = 'correlation', perEachOtherCase = False, metricsFile = 'metricsFile.h5', seed = None, PCNpath = 'data/'):

        '''Function called automatically and sets up working directory, files, and input information'''

        if seed is None:
            np.random.seed()
        else:
            np.random.seed(seed)

        self.workingDir = workingDir

        if not os.path.exists(self.workingDir):
            os.makedirs(self.workingDir)

        #if os.path.isfile('_locked') or os.path.isfile('_completed'):
        #    self.locked = True
        #else:
        #    with open('_locked', 'w'):
        #        os.utime('_locked')
        #os.remove('_locked')

        self.otherCaseDir = otherCaseDir

        if not os.path.exists(self.otherCaseDir):
            os.makedirs(self.otherCaseDir)

        self.PCNpath = PCNpath

        self.genesOfInterest = genesOfInterest
        self.knownRegulators = knownRegulators
        self.nCPUs = nCPUs
        self.panels = panels
        self.nBootstrap = nBootstrap
        self.majorMetric = majorMetric

        self.perEachOtherCase = perEachOtherCase

        self.metricsFile = metricsFile

        self.byBatchesDir = self.workingDir + 'byBatches/'
        self.bootstrapDir = self.workingDir + 'bootstrap/'
        self.dataSaveName = self.workingDir + 'data.h5'
        self.dendroDataName = self.workingDir + 'bootstrap_experiments_dendro_data.h5'

        self.bootstrapExperiments = list(range(0, nBootstrap))

        return

    standardPanels = [
        'fraction', 
        'top50', 
        'binomial', 
        'markers', 
        ]

    deprecatedPanels = [                    
        'PubMedHits', 
        'gAbove50_PanglaoMouse', 
        'gAbove50_PanglaoHuman', 
        'GOpositive', 
        'GOnegative', 
        'markerstop50', 
        'expression', 
        'closeness',
        'age', 
        'rate', 
        ]

    combinationPanels = [
        'combo3avgs',
        'combo4avgs',
        #'combo3avgs-peak',
        ]

    combinationPanelsDict = {
        'combo2avgs': ['fraction', 'binomial'],
        'combo3avgs': ['fraction', 'top50', 'binomial'],
        'combo4avgs': ['fraction', 'top50', 'binomial', 'markers'],
        'combo5-1avgs': ['fraction', 'top50', 'binomial', 'markers', 'PubMedHits'],
        'combo5-2avgs': ['fraction', 'top50', 'binomial', 'markers', 'gAbove50_PanglaoMouse'],
        'combo5-3avgs': ['fraction', 'top50', 'binomial', 'markers', 'gAbove50_PanglaoHuman'],
        'combo6-1avgs': ['fraction', 'top50', 'binomial', 'markers', 'PubMedHits', 'gAbove50_PanglaoMouse'],
        'combo6-2avgs': ['fraction', 'top50', 'binomial', 'markers', 'PubMedHits', 'gAbove50_PanglaoHuman'],
        'combo6-3avgs': ['fraction', 'top50', 'binomial', 'markers', 'gAbove50_PanglaoMouse', 'gAbove50_PanglaoHuman'],
        'combo7avgs': ['fraction', 'top50', 'binomial', 'markers', 'PubMedHits', 'gAbove50_PanglaoMouse', 'gAbove50_PanglaoHuman'],
        }

    def prepareDEG(self, dfa, dfb):

        '''Create rank dataframe with genes ranked by differential expression 
        
        Parameters:
            dfa: pandas.Dataframe
                Dataframe containing expression data for samples of interest
                Has genes as rows and batches and cells as columns 

            dfb: pandas.Dataframe
                Dataframe containing expression data for _________________
                Has genes as rows and batches and cells as columns 

        Returns:
            None 

        Usage:
            prepareDEG(dfa, dfb)
        '''
    
        print('Saving expression data', flush=True)
        dfa.to_hdf(self.dataSaveName, key='df', mode='a', complevel=4, complib='zlib')

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
        df_ranks.to_hdf(self.dataSaveName, key='df_ranks', mode='a', complevel=4, complib='zlib')
        print(df_ranks)

        return

    def preparePerBatchCase(self, **kwargs):

        ''' TEXT
        
        Parameters:
            Any parameters that function 'analyzeCase' can accept

        Returns:
            None 

        Usage: 
            an = Analysis()

            an.preparePerBatchCase()
        '''

        self.analyzeCase(pd.read_hdf(self.dataSaveName, key='df'), suffix='all', saveDir=self.byBatchesDir, toggleCalculateMajorMetric=True, toggleExportFigureData=True, toggleCalculateMeasures=True, toggleGroupBatches=False, toggleAdjustText=False, noPlot=True, **kwargs)

        return

    def _forPrepareBootstrapExperiments(self, args):

        ''' TEXT 

        Parameters: 
            saveSubDir: str
                Subdirectory for each bootstrap experiment

            df_ranks: pandas.DataFrame
                Genes ranked by differential expression 
                If None function will use rank dataframe from working directory

            df_measure: pandas.DataFrame
                Gene expression distance metric for each batch 

            df_fraction: pandas.DataFrame
                Fraction of cells expressing each gene for each batch 

            df_median_expr: pandas.DataFrame
                Median of non-zero values of each gene expression for each batch

            se_count: pandas.DataFrame
                Per batch counts 

        Returns: 
            None 

        Usage:
            self._forPrepareBootstrapExperiments((saveSubDir, df_ranks, df_measure, df_fraction, df_median_expr, se_count))

        '''

        saveSubDir, df_ranks, df_measure, df_fraction, df_median_expr, se_count = args

        np.random.seed()

        try:
            print(saveSubDir, end='\t', flush=True)
            if not os.path.exists(os.path.join(self.bootstrapDir, saveSubDir)):
                os.makedirs(os.path.join(self.bootstrapDir, saveSubDir))

            if saveSubDir == 'All':
                batches = df_measure.columns
            else:
                batches = np.random.choice(df_measure.columns, size=len(df_measure.columns), replace=True)

            np.savetxt(os.path.join(self.bootstrapDir, saveSubDir, 'batches.txt'), batches, fmt='%s')

            df_measure_temp = pd.Series(data=np.nanmedian(df_measure[batches].values.copy(), axis=1, overwrite_input=True), index=df_measure.index)
        
            df_measure_temp = df_measure_temp.unstack(0)

            df_measure_temp.to_hdf(os.path.join(self.bootstrapDir, saveSubDir, self.metricsFile), key=self.majorMetric, mode='a', complevel=4, complib='zlib')

            df_fraction_temp = df_fraction[batches]
            df_fraction_temp.columns = df_fraction_temp.columns + '_' + np.array(range(len(df_fraction_temp.columns))).astype(str)
            df_fraction_temp = df_fraction_temp.mean(axis=1)
            df_fraction_temp.to_hdf(os.path.join(self.bootstrapDir, saveSubDir, 'perGeneStats.h5'), key='df_fraction', mode='a', complevel=4, complib='zlib')

            df_median_expr_temp = df_median_expr[batches]
            df_median_expr_temp.columns = df_median_expr_temp.columns + '_' + np.array(range(len(df_median_expr_temp.columns))).astype(str)
            df_median_expr_temp = df_median_expr_temp.mean(axis=1)
            df_median_expr_temp.to_hdf(os.path.join(self.bootstrapDir, saveSubDir, 'perGeneStats.h5'), key='df_expression', mode='a', complevel=4, complib='zlib')

            se_count_temp = se_count[batches]
            se_count_temp.index = se_count_temp.index + '_' + np.array(range(len(se_count_temp.index))).astype(str)
            se_count_temp = se_count_temp.sort_values(ascending=False)
            se_count_temp.to_hdf(os.path.join(self.bootstrapDir, saveSubDir, 'perGeneStats.h5'), key='se_count', mode='a', complevel=4, complib='zlib')
        
            np.savetxt(os.path.join(self.bootstrapDir, saveSubDir, 'size.txt'), [df_fraction_temp.shape[0], se_count_temp.sum()], fmt='%i')

            if not df_ranks is None:
                df_ranks_temp = df_ranks[batches]
                df_ranks_temp.columns = df_ranks_temp.columns + '_' + np.array(range(len(df_ranks_temp.columns))).astype(str)
                df_ranks_temp = df_ranks_temp.median(axis=1).sort_values()
                df_ranks_temp.to_hdf(os.path.join(self.bootstrapDir, saveSubDir, 'perGeneStats.h5'), key='df_ranks', mode='a', complevel=4, complib='zlib')

        except Exception as exception:
            print(exception)

        return

    def prepareBootstrapExperiments(self, allDataToo = True, df_ranks = None, parallel = False):

        '''Prepare bootstrap experiments for all bootstraps by calculating gene statistics for each experiment 
        
        Parameters:
            allDataToo: boolean, Default True
                Whether to prepare experiment for all data as well

            df_ranks: pd.DataFrame, Default None
                Genes ranked by differential expression 
                If None function will use rank dataframe from working directory 

        Returns:
            None 

        Usage:
            an = Analysis()

            an.prepareBootstrapExperiments()
        '''

        try:
            if df_ranks is None:
                df_ranks=pd.read_hdf(self.dataSaveName, key='df_ranks')

            print('Reading precalculated data from %s' % self.byBatchesDir, flush=True)
            df_measure = pd.read_hdf(os.path.join(self.byBatchesDir, self.metricsFile), key=self.majorMetric)
            print('df_measure', '%1.1fMb'%(sys.getsizeof(df_measure) / 1024**2), flush=True)

            df_fraction =  pd.read_hdf(os.path.join(self.byBatchesDir, 'perGeneStats.h5'), key='df_fraction')
            df_median_expr =  pd.read_hdf(os.path.join(self.byBatchesDir, 'perGeneStats.h5'), key='df_expression')
            se_count =  pd.read_hdf(os.path.join(self.byBatchesDir, 'perGeneStats.h5'), key='se_count')

            saveSubDirs = ['Experiment %s' % (id + 1) for id in self.bootstrapExperiments]
            if allDataToo:
                saveSubDirs = ['All'] + saveSubDirs

            if parallel:
                pool = multiprocessing.Pool(processes=self.nCPUs)
                pool.map(self._forPrepareBootstrapExperiments, [(saveSubDir, df_ranks, df_measure, df_fraction, df_median_expr, se_count) for saveSubDir in saveSubDirs])
                pool.close()
                pool.join()
            else:
                for saveSubDir in saveSubDirs:
                    self._forPrepareBootstrapExperiments((saveSubDir, df_ranks, df_measure, df_fraction, df_median_expr, se_count))

            print(flush=True)

        except Exception as exception:
            print(exception)

        return

    def compareTwoCases(self, saveDir1, saveDir2, name1 = 'N1', name2='N2', saveName = 'saveName'):

        '''Compare gene measurements between two cases for each bootstrap experiment 
        
        Parameters:
            saveDir1: str
                Directory storing gene measurement dataframes for case 1

            saveDir2: str
                Directory storing gene measurement dataframes for case 2

            name1: str, Default 'N1'
                Phrase to append to keys of the resulting dataframe for case 1

            name2: str, Default 'N2'
                Phrase to append to keys of the resulting dataframe for case 2

            saveName: str, Default 'saveName'
                Name of file to save result dataframe to 

        Returns:
            None 

        Usage:
            an = Analysis()

            an.compareTwoCases(saveDir1, saveDir2, name1, name2, saveName)
        '''

        majorMetric = self.majorMetric

        try:
            df1 = pd.read_hdf(os.path.join(saveDir1, 'per-gene-measures-%s.h5' % majorMetric), key='df')
            df2 = pd.read_hdf(os.path.join(saveDir2, 'per-gene-measures-%s.h5' % majorMetric), key='df')
        except Exception as exception:
            print('ERROR reading h5 in compareTwoCases:', exception, '\nTrying reading XLSX format')

            try:
                df1 = pd.read_excel(os.path.join(saveDir1, 'per-gene-measures-%s.xlsx' % majorMetric), header=0, index_col=[0,1]).fillna('')
                df2 = pd.read_excel(os.path.join(saveDir2, 'per-gene-measures-%s.xlsx' % majorMetric), header=0, index_col=[0,1]).fillna('')
            except Exception as exception:
                print('ERROR reading XLSX:', exception)

                return

        n23_1 = len(np.intersect1d(np.unique(df1.index.get_level_values('gene').values), gEC23))
        n23_2 = len(np.intersect1d(np.unique(df2.index.get_level_values('gene').values), gEC23))

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

        if False:
            df_res_f = df_res.copy()
            df_res_f = df_res_f.loc[(df_res_f[('Intra-measures', 'EC23 count ' + name1 + ' %s' % n23_1)] >= 5) &
                                    (df_res_f[('Intra-measures', 'EC23 count ' + name2 + ' %s' % n23_2)] >= 5) &
                                    (df_res_f[('Intra-measures', 'AUC ' + name1)] >= 0.5) &
                                    (df_res_f[('Intra-measures', 'AUC ' + name2)] >= 0.5)]
            df_res_f.to_excel('%s_filtered.xlsx' % saveName, merge_cells=False)

        return

    def runPairOfExperiments(self, args):

        '''Analyze the case, compare it with comparison case, and find the conserved genes between the cases 
        
        Parameters:
            saveDir: str
                Directory for all bootstrap experiments

            saveSubDir: str 
                Subdirectory for each bootstrap experiment

            otherCaseDir: str
                Directory holding comparison data  

        Returns:
            None 

        Usage:
            self.runPairOfExperiments((saveDir, saveSubDir, otherCaseDir))
        '''

        saveDir, saveSubDir, otherCaseDir = args

        print(saveDir, saveSubDir, flush=True)

        try:
            comparisonName = os.path.join(saveDir, saveSubDir, 'comparison')

            self.analyzeCase(None, toggleAdjustText=False, noPlot=True,
                    suffix=saveSubDir, saveDir=os.path.join(saveDir, saveSubDir), printStages=False,
                    toggleCalculateMajorMetric=False, toggleExportFigureData=True, toggleCalculateMeasures=True)

            self.compareTwoCases(os.path.join(saveDir, saveSubDir, ''), otherCaseDir, name1='name1', name2='name2', saveName=comparisonName)

            additionalData = externalPanelsData.copy()
            additionalData.update({'conservedGenes': pd.read_excel(comparisonName + '.xlsx', index_col=1, header=0)['Inter-measures.T50_common_count']})

            self.analyzeCase(None, toggleAdjustText=False, dpi=300,
                    suffix=saveSubDir, saveDir=os.path.join(saveDir, saveSubDir), printStages=False,
                    toggleCalculateMajorMetric=False, toggleExportFigureData=True, toggleCalculateMeasures=True, externalPanelsData=additionalData) 

        except Exception as exception:
            print(exception)

        return

    def analyzeBootstrapExperiments(self):

        '''Analyze all bootstrap experiments and calculates dendrogram correlation data
        
        Parameters:
            None

        Returns:
            None 

        Usage:
            an = Analysis()

            an.analyzeBootstrapExperiment()
        '''

        saveSubDirs = ['All'] + ['Experiment %s' % (id + 1) for id in self.bootstrapExperiments]
        pool = multiprocessing.Pool(processes=self.nCPUs)
        pool.map(self.runPairOfExperiments, [(self.bootstrapDir, saveSubDir, os.path.join(self.otherCaseDir, 'bootstrap/', saveSubDir, '') if self.perEachOtherCase else self.otherCaseDir) for saveSubDir in saveSubDirs])
        pool.close()
        pool.join()

        dfs = []
        for id in self.bootstrapExperiments:
            saveSubDir = 'Experiment %s' % (id + 1)

            filePath = os.path.join(self.bootstrapDir, saveSubDir, 'dendrogram-heatmap-correlation-data.h5')
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

        dfs.to_hdf(self.dendroDataName, key='df', mode='a', complevel=4, complib='zlib')

        return

    def analyzeCombinationVariant(self, variant):

        ''' TEXT
        
        Parameters:
            variant: str
                Name of combination variant (e.g. 'Avg combo4avgs', 'Avg combo3avgs')

        Returns:
            pandas.DataFrame 
                TEXT

        Usage:
            an = Analysis()

            an.analyzeCombinationVariant(variant)
        '''

        def getPeaksLists(df_temp):

            '''Get list of peaks and experiments
            
            Parameters:
                df_temp: pandas.DataFrame
                    Dataframe containing dendrogram data for specific variant 

            Returns:
                list: 
                    All genes in peak for all experiments 

                list:
                    Lists of genes in peak for each experiment

                list
                    Experiments

            Usage: MARK
                getPeaksList(df_temp)
            '''

            peaksListsMerged = []
            peaksLists = []
            experiments = []
            for experiment in np.unique(df_temp.index.get_level_values('experiment')):
                se = df_temp.xs(experiment, level='experiment', axis=0).xs('species', level='species', axis=0)
                genesInHalfPeak = getGenesOfPeak(se)

                peaksListsMerged.extend(genesInHalfPeak)
                peaksLists.append(genesInHalfPeak)
                experiments.append(experiment)

            return peaksListsMerged, peaksLists, experiments

        print('Variant:', variant)

        df = pd.read_hdf(self.dendroDataName, key='df').fillna(0).set_index('gene', append=True).droplevel('order')[variant]
        
        variabilitySavePath = self.workingDir + '/variability_' + variant.replace(' ', '-') + '.xlsx'

        get_mean_std_cov_ofDataframe(df.unstack('gene').fillna(0.).T).to_excel(variabilitySavePath)

        writer = pd.ExcelWriter(self.workingDir + '%s bootstrap_in-peak_genes_SD.xlsx' % variant)

        try:
            listsMerged, lists, experiments = getPeaksLists(df.copy())
            listsSizes = [len(item) for item in lists]

            umL = np.unique(listsMerged, return_counts=True)
            se = pd.Series(index=umL[0], data=umL[1]/len(experiments)).sort_values(ascending=False)

            filePath = os.path.join(self.bootstrapDir + '/All/', 'dendrogram-heatmap-correlation-data.xlsx')
            peakGenesAll = getGenesOfPeak(pd.read_excel(filePath, index_col=0, header=0, sheet_name='Cluster index')[variant])
            df_res = pd.concat([se, se], axis=1, sort=False)
            df_res.iloc[:, 1] = np.where(np.isin(df_res.index.values, peakGenesAll), 1, np.nan)
            df_res.columns = ['Bootstrap', 'In all']

            try:
                dfv = pd.read_excel(variabilitySavePath, header=0, index_col=0).reindex(df_res.index)
                df_res = pd.concat([df_res, dfv], axis=1, sort=False)
            except Exception as exception:
                print('ERROR reading variability file:', exception)

            df_res.to_excel(writer, 'Sheet1')

            df_exp = pd.DataFrame(index=range(1000), columns=experiments)
            for i, col in enumerate(df_exp.columns):
                df_exp.iloc[:len(lists[i]), i] = lists[i]
            df_exp.to_excel(writer, 'Bootstrap lists', index=False)

        except Exception as exception:
            print('ERROR:', exception)

        writer.save()

        df_res['Bootstrap'].to_excel(self.workingDir + variant + '_variant.xlsx')

        return df_res

    def _forScramble(self, args):

        '''TEXT
        
        Parameters:
            workingDir: str 
                Directory to retrieve and save file and results to

            measures: list
                Measures (e.g: [Markers', 'Binomial -log(pvalue)', 'Top50 overlap'])

            df: pandas.DataFrame 
                TEXT 

            N: int 
                Number of iterations for _____

            maxDistance: int 
                Maximum distance away considered to be in peak

            halfWindowSize: int
                TEXT

            j: int
                Used to set randomness deterministic

        Returns:
            None

        Usage:
            self._forScramble((workingDir, measurements, df, N, maxDistance, halfWindowSize, j))
        '''

        workingDir, measures, df, N, maxDistance, halfWindowSize, j = args

        np.random.seed(j)

        allGenes = df.index.values.copy()

        listsNonmerged = []
        for i in range(N):
            np.random.shuffle(allGenes)

            data = np.zeros(len(allGenes))
            for measure in measures:
                data += movingAverageCentered(normSum1(df[measure].loc[allGenes]), halfWindowSize, looped=False)

            data = movingAverageCentered(data, halfWindowSize, looped=False)

            listsNonmerged.append(getGenesOfPeak(pd.Series(index=allGenes, data=data), maxDistance=maxDistance))

        pd.Series(listsNonmerged).to_pickle(workingDir + '%s' % j)

        return

    def scramble(self, measures, subDir = '', N = 10**4, M = 20):

        '''Randomize and prepare results for dataframe, plot count distribution, and check for variation
        
        Parameters:

            measures: list
                Measures (e.g: [Markers', 'Binomial -log(pvalue)', 'Top50 overlap'])

            subDir: str, Default ''
                Subdirectory to save dataframe to 

            N: int 
                Number of iterations for _____

            M: int 
                TEXT

        Returns:
            None

        Usage: 
            an = Analysis()

            an.scramble (measures)
        '''

        df = pd.read_excel(os.path.join(self.bootstrapDir + 'All/dendrogram-heatmap-correlation-data.xlsx'), index_col=0, header=0, sheet_name='Cluster index')

        workingDir = self.workingDir + 'random/'
        if subDir != '':
            workingDir = os.path.join(workingDir, subDir)

        if not os.path.exists(workingDir):
            os.makedirs(workingDir)

        # Run the randomization and prepare results dataframe
        if True:
            print('\nCalculating chunks', flush=True)
            pool = multiprocessing.Pool(processes=self.nCPUs)
            pool.map(self._forScramble, [(workingDir, measures, df.copy(), N, 25, 10, j) for j in range(M)])
            pool.close()
            pool.join()

            print('\nCombining chunks', flush=True)
            dfs = []
            for j in range(M):
                dfs.append(pd.read_pickle(workingDir + '%s' % j).apply(pd.Series))
                print(j, end=' ', flush=True)

            dfs = pd.concat(dfs, axis=0, sort=False).reset_index(drop=True).replace(np.nan, 'RemoveNaN')

            print('\nAligning genes', flush=True)
            df = pd.DataFrame(index=range(len(dfs)), columns=np.unique(dfs.values.flatten()), data=False, dtype=np.bool_).drop('RemoveNaN', axis=1)
            df[:] = (df.columns.values==dfs.values[..., None]).any(axis=1)
            print(df)

            print('\nRecording', flush=True)
            df.to_hdf(workingDir + 'combined_%s_aligned.h5' % M, key='df', mode='a', complevel=4, complib='zlib')

            for j in range(M):
                os.remove(workingDir + '%s' % j)

        # Save and plot counts distribution
        if True:
            df = pd.read_hdf(workingDir + 'combined_%s_aligned.h5' % M, key='df')

            se = df.sum(axis=0).sort_values(ascending=False)/df.shape[0]
            se.to_excel(workingDir + 'se_distribution.xlsx')
            se.hist(bins=250)
            plt.savefig(workingDir + 'se_distribution.png', dpi=300)
            plt.clf()

        # Check for variation of i-th quantile
        if False:
            df = pd.read_hdf(workingDir + 'combined_%s_aligned.h5' % M, key='df')

            q = 99.999
                
            res = dict()
            for i in range(1, 100):
                print(i, end=' ', flush=True)
                size = i*2*10**3
                res.update({size: np.percentile((df.sample(n=size, axis=0, replace=True).sum(axis=0).sort_values(ascending=False)/size).values, q)})

            se_per = pd.Series(res).sort_index()
            se_per.to_excel(workingDir + 'se_percentile_variation_%s.xlsx' % q)
            se_per.plot()
            plt.savefig(workingDir + 'se_percentile_variation_%s.png' % q, dpi=300)
            plt.clf()

        return

    def analyzeCase(self, df_expr, toggleCalculateMajorMetric = True, exprCutoff = 0.05, toggleExportFigureData = True, toggleCalculateMeasures = True, suffix = '', saveDir = '', toggleGroupBatches = True, dpi = 300, toggleAdjustText = True, figureSize=(8, 22), toggleAdjustFigureHeight=True, noPlot = False, halfWindowSize = 10, printStages = True, externalPanelsData = None, toggleIncludeHeatmap = True, addDeprecatedPanels = False):

        ''' Analyze, calculate, and generate plots for individual experiment
        
        Parameters:

            df_expr: pandas.Dataframe
                Take one species, one cluster (subset of clusters)

            toggleCalculateMajorMatric: boolean, Default True 
                Whether to calculate cdist of major metric 

            exprCutoff: float, Default 0.05 
                Cutoff for percent expression of input data

            toggleExportFigureData: boolean, Default True 
                Whether to export figure data 

            toggleCalculateMeasures: boolean, Default True  
                Whether to calculate measures 

            suffix: str, Default ''
                Name of experiment 

            saveDir: str, Default '' 
                Exerything is exported to this directory, should be unique for each dataset

            toggleGroupBatches: boolean, Default True
                TEXT

            dpi: int or 'figure', Default 300
                Resolution in dots per inch, if 'float' use figures dpi value 

            toggleAdjustText: boolean, Default True 
                Whether to use module to fix text overlap in figure

            figure_size: tuple, Default (8, 20)
                Width, height in inches

            toggleAdjustFigureHeight: boolean, Default True 
                Whether to adjust figure height 

            noPlot: boolean, Default False
                Whether to plot _______

            halfWindowSize: int, Default 10
                TEXT

            printStages: boolean, Default True
                Whether to print stage status to output
            
            externalPanelsData: dict, Default None 
                TEXT

            toggleIncludeHeatmap: boolean, Default True
                Whether to plot heatmap 
            
            addDeprecatedPanels: boolean, Default False
                Whether to include deprecated panels 

        Returns:
            None

        Usage: 
            self.analyzeCase(df_expr)
        '''

        stimulators, inhibitors = self.knownRegulators, []

        def calculateMajorMetricAndGeneStats(df_expr, saveDir, groupBatches, selGenes, exprCutoff):

            ''' Calculate cdist of metric (e.g. correlation)
                Calculate fraction of cells expressing each gene, and median of non-zero gene expression (per batch)

            Parameters:

                df_expr: pandas.DataFrame
                    Take one species, one cluster (subset of clusters)

                saveDir: str
                    Exerything is exported to this directory, should be unique for each dataset

                groupBatches: boolean
                    Whether to take median across batches

                selGenes: list 
                    List of receptors, or transcription factors

                exprCutoff: float
                    Cutoff for percent expression of input data

            Returns:
                None 

            Usage:
                calculateMajorMetricAndGeneStats(df_expr, saveDir, groupBatches, selGenes, exprCutoff)
            '''

            print('Received expression data of shape:', df_expr.shape, flush=True)
            np.savetxt(os.path.join(saveDir, 'size.txt'), np.array(df_expr.shape), fmt='%i')

            # For each batch calculate gene expression distance metric
            print('Calculating distance metric', flush=True)
            df_measure = get_df_distance(df_expr, metric=self.majorMetric, genes=selGenes, analyzeBy='batch', minSize=10, groupBatches=groupBatches, cutoff=exprCutoff, nCPUs=self.nCPUs)

            print('Recording major metric (shape: %s, %s) to h5' % df_measure.shape, flush=True)
            df_measure.to_hdf(os.path.join(saveDir, self.metricsFile), key=self.majorMetric, mode='a', complevel=4, complib='zlib')

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

            ''' Builds and plots dendrogram, heatmap, and bargraphs.

            Parameters:
                df: pandas.DataFrame 
                    All calculated measurement 

                metric: str, Default 'euclidean'
                    Name of metric used to build dendrogram and identify clusters in it
                    Metric has to be of type "Euclidean" to use linkage method "Ward"
                    With any other metric (e.g. correlation distance) use linkage method "average" etc.
                    Metric 'euclidean_missing' used commonly-non-missing points only

                linkageMethod: str, Default 'ward'
                    Linkage algorithm to use

                n_clusters: int, Default 10 
                    Specific number of clusters to find 

                adjustText: str, Default toggleAdjustText
                    Whether to use module to fix text overlap in figure

            Returns:
                dict:
                    Figure data for export 

            Usage: 
                makeCombinationPlot(df)
            '''
             
            nonlocal figureSize, self

            if self.panels is None:
                self.panels = self.combinationPanels + self.standardPanels

                if addDeprecatedPanels:
                    self.panels += deprecatedPanels

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
                        adjust_text(texts, va='top', ha='center', autoalign='x', lim=400, only_move={'text':'x'}, force_text=(0.2, 0.5))

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

                    clb = fig.colorbar(im, ax=ax, fraction=0.4, label='Eucl. dist. of gene expr.\n %s dist.' % self.majorMetric)
                    clb.ax.tick_params(labelsize=fontsize)

                return

            def addBar(fig, dataArgs, panel, coords, halfWindowSize = halfWindowSize, noAverage = False):

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
                showBoundaries = True

                # Draw randomized bar panels
                if False:
                    showBoundaries = False
                    np.random.seed(0)
                    rorder = np.arange(len(allGenes))
                    np.random.shuffle(rorder)
                    allGenes = allGenes[rorder]
                    tickLabelsColors = tickLabelsColors[rorder]

                ax = fig.add_axes(coords, frame_on=True)
                ax.set_xlim([min(clusterBoundaries), max(clusterBoundaries)])

                if panel == 'fraction':
                    ylabel='Fraction'
                    try:
                        data = pd.read_hdf(os.path.join(saveDir, 'perGeneStats.h5'), key='df_fraction').reindex(allGenes)
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
                    ylabel='PubMedHits'
                    try:
                        data = pd.Series(externalPanelsData['pubMed angiogenesis hits']).reindex(allGenes).values
                        #data = np.log(data)
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

                        try:
                            data = data[rorder]
                        except:
                            pass

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
                        data = -np.log(binomialEnrichmentProbability('PCN.txt', enriched_genes=diffExpressedGenes, target_genes=selGenes, PCNpath=self.PCNpath)['Binomial_Prob'].reindex(allGenes).values)
                    except:
                        data = np.zeros(len(allGenes))

                        try:
                            data = externalPanelsData['diffExpressedGenes']
                            if not type(data) is pd.Series:
                                data = data.median(axis=1)

                            data = data[data > 0].sort_values()[:1000]
                            diffExpressedGenes = data.index

                            data = -np.log(binomialEnrichmentProbability('PCN.txt', enriched_genes=diffExpressedGenes, target_genes=selGenes, PCNpath=self.PCNpath)['Binomial_Prob'].reindex(allGenes).values)
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

                elif panel[-len('-peak'):] == '-peak':
                    ylabel = panel
                    try:
                        targetPanel = 'Avg ' + panel[:-len('-peak')]
                        data = pd.Series(index=allGenes, data=np.nan_to_num(panelsData[panelsDataNames[targetPanel]]))
                        data /= data.max()
                        peaks = getGenesOfPeak(data)
                        data[:] = 0.
                        data[peaks] = 1.

                        noAverage = True
                    except:
                        data = np.zeros(len(allGenes))

                elif panel[:len('combo')] == 'combo':
                    ylabel = panel
                    try:
                        def norm1(s):

                            w = np.nansum(panelsData[panelsDataNames[s]])
                            if w != w or w == 0.:
                                w = 1.

                            return np.nan_to_num(panelsData[panelsDataNames[s]]) / w

                        data = np.zeros(len(allGenes))
                        for subPanel in self.combinationPanelsDict[panel]:
                            if panel[-len('avgs'):] == 'avgs':
                                data += movingAverageCentered(norm1(subPanel), halfWindowSize)
                            else:
                                data += norm1(subPanel)
                    except:
                        data = np.zeros(len(allGenes))

                ax.bar(range(len(clusters)), data, width=ax.get_xlim()[1]/len(clusters), color=tickLabelsColors)

                data_avg = movingAverageCentered(np.nan_to_num(data), halfWindowSize)

                if not noAverage:
                    ax.plot(range(len(clusters)), data_avg, linewidth=1.0, color='coral', alpha=1.0)
                    ax.text(0.999, 0.95, 'window = %s' % (2*halfWindowSize + 1), color='coral', ha='right', va='top', transform=ax.transAxes, fontsize=3)

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

                if showBoundaries:
                    ylim = ax.get_ylim()
                    for i in range(1, len(np.unique(clusters))):
                        ax.plot([clusterBoundaries[i] - 0.5]*2, [ylim[0], ylim[1]], '--', lw=0.5, color='k', clip_on=False)

                ax.text(-0.01, 0.5, ylabel, fontsize=4, rotation=0, va='center', ha='right', transform=ax.transAxes)

                return ylabel, data, data_avg

            if metric == 'euclidean_missing':
                metric = metric_euclidean_missing

            mmin, mmax = np.nanmin(np.nanmin(df.values)), np.nanmax(np.nanmax(df.values))

            if self.majorMetric == 'correlation':
                missingFillValue = 1.0
            elif self.majorMetric == 'cosine':
                missingFillValue = 1.0
            elif self.majorMetric == 'euclidean':
                missingFillValue = mmax
            else:
                missingFillValue = mmax

            if printStages:
                print('\tCalculating %s metric of %s . . .' % (metric, self.majorMetric), end='\t', flush=True)
            M = pdist(df.fillna(missingFillValue).values.T, metric=metric)
            if printStages:
                print('Done', flush=True)

            nPanels = len(self.panels)
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
                for ipanel, panel in enumerate(reversed(self.panels)):
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

                number = np.loadtxt(os.path.join(saveDir, 'size.txt'), dtype=int)
                fig.suptitle('Ordering by single cell co-expression\nData: %s\n(%s cells, %s x %s genes)' % (suffix, number[1], number[0], df.shape[1]), fontsize=8, backgroundcolor='white')

                if printStages:
                    print('\tSaving image . . .', end='\t', flush=True)
                fig.savefig(os.path.join(saveDir, '%s dendrogram-heatmap-%s.png' % (suffix, self.majorMetric)), dpi=dpi)
                if printStages:
                    print('Done', flush=True)

            plt.close(fig)

            return dataArgs

        def exportFigureData(dataArgs, saveXLSX = True, saveHDF = True):

            ''' Function to assist in exporting figure data to excel and/or hdf file 

            Parameters:

                dataArgs: dict
                    Figure data for export 

                saveXLSX: boolean, Default True 
                    Whether to export data as an excel file

                saveHDF: boolean, Default True 
                    Whether to export data as an hdf file

            Returns:
                None 

            Usage:
                exportFigureData(dataArgs)
            '''

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
                writer = pd.ExcelWriter(os.path.join(saveDir, 'dendrogram-heatmap-%s-data.xlsx' % self.majorMetric))
                df_C.to_excel(writer, 'Cluster index')
                df_M.to_excel(writer, 'Expression distance measure')
                writer.save()

            if saveHDF:
                df_C.to_hdf(os.path.join(saveDir, 'dendrogram-heatmap-%s-data.h5' % self.majorMetric), key='df_C', mode='a', complevel=4, complib='zlib')
                df_M.to_hdf(os.path.join(saveDir, 'dendrogram-heatmap-%s-data.h5' % self.majorMetric), key='df_M', mode='a', complevel=4, complib='zlib')

            return

        if not os.path.exists(saveDir):
            os.makedirs(saveDir)

        selGenes = np.unique(list(self.genesOfInterest))

        if toggleCalculateMajorMetric and (not df_expr is None):
            calculateMajorMetricAndGeneStats(df_expr, saveDir, toggleGroupBatches, selGenes, exprCutoff)

        # Load and prepare df_measure
        if True:
            try:
                df_measure = pd.read_hdf(os.path.join(saveDir, self.metricsFile), key=self.majorMetric)

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
            df_Em = pd.read_hdf(os.path.join(saveDir, 'dendrogram-heatmap-%s-data.h5' % self.majorMetric), key='df_M')
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

                try:
                    AUC = roc_auc_score(~np.isin(se.index.values, selGenes23), se.values) if len(se) >= 10 else np.nan
                except:
                    AUC = np.nan

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

            df.to_excel(os.path.join(saveDir, 'per-gene-measures-%s.xlsx' % self.majorMetric), merge_cells=False)
            df.to_hdf(os.path.join(saveDir, 'per-gene-measures-%s.h5' % self.majorMetric), key='df', mode='a', complevel=4, complib='zlib')
            if printStages:
                print('Done', flush=True)

        return

    def reanalyzeMain(self, **kwargs):

        ''' Reanalyze all-data case 

        Parameters: 
            Any parameters that function 'analyzeCase' can accept

        Returns:
            None 

        Usage:
            an = Analyze()

            an.reanalyzeMain()
        '''

        try:
            print('Re-analyzing all-data case', flush=True)
            
            self.analyzeCase(None, toggleAdjustText=True, dpi=300, suffix='All', saveDir=os.path.join(self.bootstrapDir, 'All/'), printStages=True, toggleCalculateMajorMetric=False, toggleExportFigureData=True, toggleCalculateMeasures=False, externalPanelsData=dict(externalPanelsData, conservedGenes=pd.read_excel(os.path.join(self.bootstrapDir, 'All', 'comparison.xlsx'), index_col=1, header=0)['Inter-measures.T50_common_count']), **kwargs)

            shutil.copyfile(os.path.join(self.bootstrapDir, 'All', 'All dendrogram-heatmap-correlation.png'), os.path.join(self.workingDir, 'results.png'))

        except Exception as exception:
            print(exception)

        return

    def analyzeAllPeaksOfCombinationVariant(self, variant, nG = 8, nE = 30, fcutoff = 0.5, width = 50):

        ''' Reanalyze all-data case 

        Parameters: 
            variant: str
                Name of combination variant (e.g. 'Avg combo4avgs', 'Avg combo3avgs')

            nG: int, Default 8
                Threshold to apply when forming flat clusters MARK

            nE: int, Default 30
                Threshold to apply when forming flat clusters MARK

            fcutoff: float, Default 0.5
                Fraction cutoff MARK

            width: int, Default 50
                Width of peak 

        Returns:
            None 

        Usage:
            an = Analyze()

            an.analyzeAllPeaksOfCombinationVariant('Avg combo3avgs', nG=8, nE=30, fcutoff=0.5, width=50)
            
            an.analyzeAllPeaksOfCombinationVariant('Avg combo4avgs', nG=8, nE=30, fcutoff=0.5, width=50)
        '''

        def getPeaksLists(df_temp):

            ''' Get list of peaks and genes

            Parameters:
                df_temp: pandas.DataFrame
                    Dendrogram data

            Returns:
                list:
                    Lists of genes in peak for each experiment

                ndarray:
                    Sorted and unique genes 

            Usage:
                getPeaksList(df_temp)
            '''

            peaksLists = {}
            genes = []
            for experiment in np.unique(df_temp.index.get_level_values('experiment')):
                se = df_temp.xs(experiment, level='experiment', axis=0)
                peaksLists.update({(experiment, i): getGenesOfPeak(se, peak=peak, maxDistance=int(width/2)) for i, peak in enumerate(getPeaks(se, distance=width))})
                genes.extend(se.index.values.tolist())

            return peaksLists, np.unique(genes)

        print('Variant:', variant)

        df = pd.read_hdf(self.dendroDataName, key='df').fillna(0).set_index('gene', append=True).droplevel('order')[variant]

        listsDict, allgenes = getPeaksLists(df.xs('species', level='species', axis=0).copy())

        se_peakAssignments = pd.Series(listsDict).apply(pd.Series)
        df_m = pd.DataFrame(index=allgenes, data=0, columns=se_peakAssignments.index)
        for peak in df_m.columns:
            df_m.loc[se_peakAssignments.loc[peak].dropna().values, peak] = 1

        df_m = df_m.loc[df_m.sum(axis=1) > 0]

        Z1 = hierarchy.linkage(df_m, 'ward')
        Z2 = hierarchy.linkage(df_m.T, 'ward')

        O1 = hierarchy.dendrogram(Z1, no_plot=True, get_leaves=True)['leaves']
        O2 = hierarchy.dendrogram(Z2, no_plot=True, get_leaves=True)['leaves']

        df_m = df_m.iloc[O1]
        df_m = df_m.T.iloc[O2].T

        clusters1 = hierarchy.fcluster(Z1, t=nG, criterion='maxclust')[O1] - 1
        clusters2 = hierarchy.fcluster(Z2, t=nE, criterion='maxclust')[O2] - 1

        u1 = np.unique(clusters1, return_counts=True)[1]
        u2 = np.unique(clusters2, return_counts=True)[1]

        df_m.index = pd.MultiIndex.from_arrays([df_m.index, clusters1], names=['gene', 'cluster'])
        df_m.columns = pd.MultiIndex.from_arrays(['E' + df_m.columns.get_level_values(0).str.split('Experiment ', expand=True).get_level_values(-1) + '.' + df_m.columns.get_level_values(1).astype(str), clusters2], names=['peak', 'cluster'])

        ndict = dict()
        res = dict()
        for ci in df_m.index.levels[-1]:
            for cj in df_m.columns.levels[-1]:
                df_temp = df_m.xs(ci, level='cluster', axis=0).xs(cj, level='cluster', axis=1)
                m = df_temp.values.mean()
                if m >= fcutoff:
                    #print(ci, cj, df_temp.shape, '\t', np.round(m, 2), '\t', cleanListString(sorted(df_temp.index.values.tolist())))
                    res[(ci, cj)] = df_temp.shape[1]
                    ndict[ci] = cleanListString(sorted(df_temp.index.values.tolist()))

        se = pd.Series(res).groupby(level=0).sum()
        se.name = 'frequency'

        df = se.to_frame()
        df['genes'] = pd.Series(se.index).replace(ndict).values
        df['frequency'][:] = 0

        for i, group in enumerate(df['genes'].values):
            group = group.split(', ')
            fractions = df_m.droplevel('cluster').loc[group].sum(axis=0) / len(group)
            df['frequency'].iloc[i] = len(fractions[fractions >= fcutoff]) / df_m.shape[1]

        df = df.sort_values(by='frequency', ascending=False)
        print(df)
        
        df.to_excel(self.workingDir + 'All peaks %s. nG%s-nE%s.xlsx' % (variant, nG, nE), index=False)

        return
