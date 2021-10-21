''' Module holding class that implements the analysis pipeline 
'''

from .commonFunctions import *

class Analysis():

    '''Class of analysis and visualization functions for DECNEO
    
    Parameters:
        workingDir: str, Default ''
            Directory to retrieve and save files and results to
        
        otherCaseDir: str, Default ''
            Directory holding comparison (other species) data
        
        genesOfInterest: list, Default None
            Particular genes to analyze, e.g. receptors

        knownRegulators: list, Default None
            Known marker genes 

        nCPUs: int, Default 1
            Number of CPUs to use for multiprocessing, recommended 10-20

        panels: list, Default None
            Particular measurements to include in the analysis

        nBootstrap: int, Default 100
            Number of bootstrap experiments to perform

        majorMetric: str, Default 'correlation'
            Metric name (e.g. 'correlation', 'cosine', 'euclidean', 'spearman')

        methodForDEG: str, Default 'ttest'
            Possible options: {'ttest', 'mannwhitneyu'}

        perEachOtherCase: boolean, Default False 
            Whether to perform comparisons of bootstrap experiments with other bootstrap experiments or with a single case

        metricsFile: str, 'metricsFile.h5' 
            Name of file where gene expression distance data is saved for specified metric

        seed: int, None 
            Used to set randomness deterministic

        PCNpath: str, Default 'data/'
            Path to PCN file
       
    '''

    def __init__(self, workingDir = '', otherCaseDir = '', genesOfInterest = None, knownRegulators = None, nCPUs = 1, panels = None, nBootstrap = 100, majorMetric = 'correlation', perEachOtherCase = False, metricsFile = 'metricsFile.h5', seed = None, PCNpath = 'data/', minBatches = 5, pseudoBatches = 10, dendrogramMetric = 'euclidean', dendrogramLinkageMethod = 'ward', methodForDEG = 'ttest'):

        '''Function called automatically and sets up working directory, files, and input information'''

        if seed is None:
            np.random.seed()
        else:
            np.random.seed(seed)

        self.workingDir = workingDir

        if not os.path.exists(self.workingDir):
            os.makedirs(self.workingDir)

        ## Locking mechanism is not implemented
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
        self.dendrogramMetric = dendrogramMetric
        self.dendrogramLinkageMethod = dendrogramLinkageMethod
        self.methodForDEG = methodForDEG

        self.minBatches = minBatches
        self.pseudoBatches = pseudoBatches

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
        'binomial', 
        'top50', 
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
        #'combo3avgs-peak', # Use this syntax to include panels with genes of peak
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

    def prepareDEG(self, dfa, dfb, pvalueLimit = 0.001):

        '''Save gene expression data of cell type of interest.
        Create rank dataframe (df_ranks) with genes ranked by differential expression 
        
        Parameters:
            dfa: pandas.Dataframe
                Dataframe containing expression data for cell type of interest
                Has genes as rows and (batches, cells) as columns 

            dfb: pandas.Dataframe
                Dataframe containing expression data for cells of type other than cell type of interest
                Has genes as rows and (batches, cells) as columns 

            pvalueLimit: float, Default 0.001
                Maximum possible p-value to include

        Returns:
            None 

        Usage:
            prepareDEG(dfa, dfb)
        '''
    
        dfb = dfb.reindex(dfa.index).fillna(0.)

        nOriginalBatches = len(np.unique(dfa.columns.get_level_values('batch').values))
        if nOriginalBatches < self.minBatches:
            dfa.columns = pd.MultiIndex.from_arrays([np.random.permutation(np.hstack([np.array([str(i)]*len(v), dtype=str) for i, v in enumerate(np.array_split(dfa.columns.get_level_values('batch').values, self.pseudoBatches))])), dfa.columns.get_level_values('cell')], names=['batch', 'cell'])
            dfb.columns = pd.MultiIndex.from_arrays([np.random.permutation(np.hstack([np.array([str(i)]*len(v), dtype=str) for i, v in enumerate(np.array_split(dfb.columns.get_level_values('batch').values, self.pseudoBatches))])), dfb.columns.get_level_values('cell')], names=['batch', 'cell'])

            print('Original bathces:', nOriginalBatches, '\tGenerated %s pseudobatches:' % self.pseudoBatches, dfa.shape, dfb.shape)

        print('Saving expression data', flush=True)
        dfa.to_hdf(self.dataSaveName, key='df', mode='a', complevel=4, complib='zlib')

        genes = []
        batches = []
        for batch in np.unique(dfa.columns.get_level_values('batch').values):

            try:
                df_temp_a = dfa.xs(batch, level='batch', axis=1, drop_level=False)
                df_temp_b = dfb.xs(batch, level='batch', axis=1, drop_level=False)

                if self.methodForDEG == 'ttest':
                    ttest = scipy.stats.ttest_ind(df_temp_a.values, df_temp_b.values, axis=1)
                    df_test = pd.DataFrame(index=dfa.index, columns=['statistic', 'pvalue'])
                    df_test['statistic'] = ttest[0]
                    df_test['pvalue'] = ttest[1]                            
                    df_test = df_test.sort_values('statistic', ascending=False).dropna()

                elif self.methodForDEG == 'mannwhitneyu':
                    df_temp_a = df_temp_a.apply(np.array, axis=1)
                    df_temp_b = df_temp_b.apply(np.array, axis=1)
                    df_temp = pd.concat([df_temp_a, df_temp_b], axis=1, sort=False)   
                    df_temp = df_temp.loc[(df_temp[0].apply(np.sum) + df_temp[1].apply(np.sum)) > 0]
                    df_test = df_temp.apply(lambda v: pd.Series(np.array(scipy.stats.mannwhitneyu(v[0], v[1]))), axis=1)
                    df_test.columns = ['statistic', 'pvalue']
                    df_test = df_test.sort_values('pvalue', ascending=True).dropna()
                
                genes.append(df_test.loc[df_test['pvalue'] <= pvalueLimit]['statistic'].index.values)
                batches.append(batch)

            except Exception as exception:
                print(exception)

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

        ''' Process gene expression data to generate per-batch distance measure and save to file. No plots are generated
        
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

        ''' Function used internally in multiprocessing 

        Parameters: 
            saveSubDir: str
                Subdirectory for each bootstrap experiment

            df_ranks: pandas.DataFrame
                Genes ranked by differential expression. This is a legacy parameter.
                If None function will use rank dataframe from working directory

            df_measure: pandas.DataFrame
                Gene expression per-batch distance

            df_fraction: pandas.DataFrame
                Fraction of cells expressing each gene for each batch 

            df_median_expr: pandas.DataFrame
                Median of non-zero values of each gene expression for each batch

            se_count: pandas.DataFrame
                Per batch counts 

        Returns: 
            None 

        Usage:
            For internal use only

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
                batches = np.loadtxt(os.path.join(self.bootstrapDir, saveSubDir, 'batches.txt'), dtype=str, delimiter='\t')

                df_ranks_temp = df_ranks[df_ranks.columns.intersection(batches)]
                df_ranks_temp.columns = df_ranks_temp.columns + '_' + np.array(range(len(df_ranks_temp.columns))).astype(str)
                df_ranks_temp = df_ranks_temp.median(axis=1).sort_values()
                df_ranks_temp.to_hdf(os.path.join(self.bootstrapDir, saveSubDir, 'perGeneStats.h5'), key='df_ranks', mode='a', complevel=4, complib='zlib')

        except Exception as exception:
            print(exception)

        return

    def prepareBootstrapExperiments(self, allDataToo = True, df_ranks = None, parallel = False):

        '''Prepare bootstrap experiments data and calculating gene statistics for each experiment 
        
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
                    #self._forPrepareBootstrapExperiments((saveSubDir, df_ranks, None, None, None, None))

            print(flush=True)

        except Exception as exception:
            print(exception)

        return

    def compareTwoCases(self, saveDir1, saveDir2, name1 = 'N1', name2='N2', saveName = 'saveName'):

        '''Compare gene measurements between two cases for each bootstrap experiment 
        
        Parameters:
            saveDir1: str
                Directory storing gene measurement data for case 1

            saveDir2: str
                Directory storing gene measurement data for case 2

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

        n23_1 = len(np.intersect1d(np.unique(df1.index.get_level_values('gene').values), gEC22))
        n23_2 = len(np.intersect1d(np.unique(df2.index.get_level_values('gene').values), gEC22))

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

        '''Analyze the case, compare it with comparison case, find the conserved genes between the cases, analyze case again
        
        Parameters:
            saveDir: str
                Directory with all bootstrap experiments

            saveSubDir: str 
                Subdirectory for a bootstrap experiment

            otherCaseDir: str
                Directory holding comparison data  

        Returns:
            None 

        Usage:
            For internal use only
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

        '''Analyze all bootstrap experiments
        
        Parameters:
            None

        Returns:
            None 

        Usage:
            an = Analysis()

            an.analyzeBootstrapExperiments()
        '''

        saveSubDirs = ['All'] + ['Experiment %s' % (id + 1) for id in self.bootstrapExperiments]
        pool = multiprocessing.Pool(processes=self.nCPUs)
        pool.map(self.runPairOfExperiments, [(self.bootstrapDir, saveSubDir, os.path.join(self.otherCaseDir, 'bootstrap/', saveSubDir, '') if self.perEachOtherCase else self.otherCaseDir) for saveSubDir in saveSubDirs])
        pool.close()
        pool.join()

        dfs = []
        for id in self.bootstrapExperiments:
            saveSubDir = 'Experiment %s' % (id + 1)

            filePath = os.path.join(self.bootstrapDir, saveSubDir, 'dendrogram-heatmap-%s-data.h5' % self.majorMetric)
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

        try:
            dfs = pd.concat(dfs, axis=0, sort=False)
            print(dfs)

            dfs.to_hdf(self.dendroDataName, key='df', mode='a', complevel=4, complib='zlib')
        except Exception as exception:
            print(exception)
            pass

        return

    def analyzeCombinationVariant(self, variant):

        ''' Analyze a combination of measures (same as in panels)
        
        Parameters:
            variant: str
                Name of combination variant (e.g. 'Avg combo4avgs', 'Avg combo3avgs')

        Returns:
            pandas.DataFrame 
                Analysis result

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

            Usage:
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

            filePath = os.path.join(self.bootstrapDir + '/All/', 'dendrogram-heatmap-%s-data.xlsx' % self.majorMetric)
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

        '''Function used internally in multiprocessing
        
        Parameters:
            workingDir: str 
                Working directory to retrieve and save file and results to

            measures: list
                Measures (e.g: [Markers', 'Binomial -log(pvalue)', 'Top50 overlap'])

            df: pandas.DataFrame 
                Ordered genes data

            N: int 
                Size of a chunk

            maxDistance: int 
                Maximum distance away considered to be in peak

            halfWindowSize: int
                Moving average half-window size

            j: int
                Identifier of a chunk

        Returns:
            None

        Usage:
            For internal use only
        '''

        workingDir, measures, df, N, maxDistance, halfWindowSize, j, getMax = args

        np.random.seed()

        allGenes = df.index.values.copy()

        listsNonmerged = []
        for i in range(N):
            np.random.shuffle(allGenes)

            data = np.zeros(len(allGenes))
            for measure in measures:
                data += movingAverageCentered(normSum1(df[measure].loc[allGenes]), halfWindowSize, looped=False)

            data = movingAverageCentered(data, halfWindowSize, looped=False)

            listsNonmerged.append(getGenesOfPeak(pd.Series(index=allGenes, data=data), maxDistance=maxDistance))

        se = pd.Series(listsNonmerged)

        if getMax:
            dfs = se.apply(pd.Series).reset_index(drop=True).replace(np.nan, 'RemoveNaN')
            df = pd.DataFrame(index=range(len(dfs)), columns=np.unique(dfs.values.flatten()), data=False, dtype=np.bool_)
            try:
                df = df.drop('RemoveNaN', axis=1)
            except:
                pass

            df[:] = (df.columns.values==dfs.values[..., None]).any(axis=1)

            return df.mean(axis=0).max()

        else:
            se.to_pickle(workingDir + '%s' % j)

        return

    def scramble(self, measures, subDir = '', case='All', N = 10**4, M = 20, getMax = False, maxSuff = ''):

        '''Run control analysis for the dendrogram order
        
        Parameters:
            measures: list
                Measures (e.g: [Markers', 'Binomial -log(pvalue)', 'Top50 overlap'])

            subDir: str, Default ''
                Subdirectory to save dataframe to 

            N: int 
                Chunk size

            M: int 
                Number of chunks

        Returns:
            None

        Usage: 
            an = Analysis()

            an.scramble (measures)
        '''

        df = pd.read_excel(os.path.join(self.bootstrapDir + '%s/dendrogram-heatmap-%s-data.xlsx' % (case, self.majorMetric)), index_col=0, header=0, sheet_name='Cluster index')

        workingDir = self.workingDir + 'random/'
        if subDir != '':
            workingDir = os.path.join(workingDir, subDir)

        if not os.path.exists(workingDir):
            os.makedirs(workingDir)

        # Run the randomization and prepare results dataframe
        if True:
            print('\nCalculating chunks', flush=True)
            pool = multiprocessing.Pool(processes=self.nCPUs)
            result = pool.map(self._forScramble, [(workingDir, measures, df.copy(), N, 25, 10, j, getMax) for j in range(M)])
            pool.close()
            pool.join()

            if getMax:
                pd.Series(result).to_hdf(workingDir + 'max%s%s.h5' % (N, maxSuff), key='df', mode='a', complevel=4, complib='zlib')

                return

            print('\nCombining chunks', flush=True)
            dfs = []
            for j in range(M):
                dfs.append(pd.read_pickle(workingDir + '%s' % j).apply(pd.Series))
                print(j, end=' ', flush=True)

            dfs = pd.concat(dfs, axis=0, sort=False).reset_index(drop=True).replace(np.nan, 'RemoveNaN')

            print('\nAligning genes', flush=True)
            df = pd.DataFrame(index=range(len(dfs)), columns=np.unique(dfs.values.flatten()), data=False, dtype=np.bool_)
            try:
                df = df.drop('RemoveNaN', axis=1)
            except:
                pass

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

        #os.remove(workingDir + 'combined_%s_aligned.h5' % M)

        return

    def analyzeCase(self, df_expr, toggleCalculateMajorMetric = True, exprCutoff = 0.05, toggleExportFigureData = True, toggleCalculateMeasures = True, suffix = '', saveDir = '', toggleGroupBatches = True, dpi = 300, toggleAdjustText = True, markersLabelsRepelForce = 1.5, figureSize=(8, 22), toggleAdjustFigureHeight=True, noPlot = False, halfWindowSize = 10, printStages = True, externalPanelsData = None, toggleIncludeHeatmap = True, addDeprecatedPanels = False, includeClusterNumber = True, togglePublicationFigure = False):

        '''Analyze, calculate, and generate plots for individual experiment
        
        Parameters:
            df_expr: pandas.Dataframe
                Gene expression data

            toggleCalculateMajorMatric: boolean, Default True 
                Whether to calculate cdist of major metric. This is a legacy parameter 

            exprCutoff: float, Default 0.05 
                Cutoff for percent expression in a batch of input data

            toggleExportFigureData: boolean, Default True 
                Whether to export figure data 

            toggleCalculateMeasures: boolean, Default True  
                Whether to calculate measures 

            suffix: str, Default ''
                Name of experiment 

            saveDir: str, Default '' 
                Exerything is exported to this directory, should be unique for each dataset

            toggleGroupBatches: boolean, Default True
                Whether to group batches or save per-batch distance measure

            dpi: int or 'figure', Default 300
                Resolution in dots per inch, if 'float' use figures dpi value 

            toggleAdjustText: boolean, Default True 
                Whether to use (external) module to minimize text overlap in figure

            figure_size: tuple, Default (8, 20)
                Width, height in inches

            toggleAdjustFigureHeight: boolean, Default True 
                Whether to adjust figure height 

            noPlot: boolean, Default False
                Whether to generate plot

            halfWindowSize: int, Default 10
                Moving average half-window size

            printStages: boolean, Default True
                Whether to print stage status to output
            
            externalPanelsData: dict, Default None 
                Dictionary containing additional panels data

            toggleIncludeHeatmap: boolean, Default True
                Whether to include heatmap in figure
            
            addDeprecatedPanels: boolean, Default False
                Whether to include deprecated panels 

        Returns:
            None

        Usage: 
            self.analyzeCase(df_expr)
        '''

        stimulators, inhibitors = self.knownRegulators, []

        #if togglePublicationFigure:
        #    toggleExportFigureData = True

        def calculateMajorMetricAndGeneStats(df_expr, saveDir, groupBatches, selGenes, exprCutoff):

            '''Calculate cdist of metric (e.g. correlation)
                Calculate fraction of cells expressing each gene, and median of non-zero gene expression (per batch)

            Parameters:
                df_expr: pandas.DataFrame
                    Gene expression of one species, one cluster (subset of clusters)

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

        def makeCombinationPlot(df, n_clusters = 10, adjustText = toggleAdjustText):

            '''Builds and plots dendrogram, heatmap, and bargraphs.

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
             
            nonlocal figureSize, self, togglePublicationFigure, markersLabelsRepelForce, includeClusterNumber

            metric = self.dendrogramMetric
            linkageMethod = self.dendrogramLinkageMethod
            
            if metric == 'euclidean_missing':
                metric = metric_euclidean_missing

            if self.panels is None:
                self.panels = self.standardPanels

                if addDeprecatedPanels:
                    self.panels += deprecatedPanels

                self.panels += self.combinationPanels

            def addDendro(fig, dataGenes, M, coords, metric = metric, linkageMethod = 'ward', linewidth = 0.25, adjustText = adjustText, fontsize = 6):

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
                        adjustTexts1D(texts, fig, ax)

                        #adjust_text(texts, va='top', ha='center', autoalign='x', lim=400, only_move={'text':'x'}, force_text=(markersLabelsRepelForce, 0.5))

                    v = 0.05 * ax.get_ylim()[1]
                    for text, opos in zip(texts, origPos):
                        text._y = -v
                        ax.plot([text._x, opos], [text._y, 0.], color=text._color, lw=0.5, clip_on=False)

                if True:
                    clusters = scipy.cluster.hierarchy.fcluster(Z, t=n_clusters, criterion='maxclust')[D['leaves']] - 1
                    clusterBoundaries = (np.where(clusters - np.roll(clusters, 1) != 0)[0]/ len(D['leaves'])) * ax.get_xlim()[1]
                    clusterBoundaries = np.append(clusterBoundaries, ax.get_xlim()[1])

                    clusterCenters = clusterBoundaries[:-1] + ((clusterBoundaries - np.roll(clusterBoundaries, 1))/2.)[1:]
                    vposition = (Z[-n_clusters,2] + Z[-n_clusters+1,2]) / 5

                    if includeClusterNumber:
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

            def addHeatmap(fig, dataArgs, coords, adjustText = adjustText, fontsize = 6):

                plottingMajorMetricOfSelected = False
                
                if plottingMajorMetricOfSelected:
                    M = dataArgs['majorMetricOfSelected']
                else:
                    M = dataArgs['M'] 

                order = dataArgs['order'] 
                genes = dataArgs['genes'] 
                locations = dataArgs['locations'] 
                tickLabelsColors =  dataArgs['tickLabelsColors'] 
                tickslabels =  dataArgs['xtickslabels'] 
                clusters =  dataArgs['clusters'] 
                clusterBoundaries =  dataArgs['clusterBoundaries']
                clusterCenters =  dataArgs['clusterCenters']

                ax = fig.add_axes(coords, frame_on=True)

                masked_M = np.ma.array(M, mask=np.isnan(M))

                if plottingMajorMetricOfSelected:
                    cmap = copy.copy(plt.cm.bwr)
                    cmap.set_bad('grey')
                    vmin, vmax = -1, 1
                else:
                    cmap = copy.copy(plt.cm.bwr) # Greens_r
                    cmap.set_bad('red')
                    vmin, vmax = None, None

                im = ax.imshow(masked_M, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax, interpolation='None', extent=(-0.5, M.shape[0] - 0.5, M.shape[1] - 0.5, -0.5))

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
                        adjustTexts1D(texts, fig, ax)

                        #adjust_text(texts, va='top', ha='center', autoalign='x', lim=400, only_move={'text':'x'}, force_text=(1.05*markersLabelsRepelForce, 0.5))

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
                        adjustTexts1D([ele for ele in reversed(texts)], fig, ax)

                        #adjust_text(texts, va='center', ha='right', autoalign='y', lim=400, only_move={'text':'y'}, force_text=(1.05*markersLabelsRepelForce, 0.5))

                    v = -0.01 * ax.get_xlim()[1]
                    for text, opos in zip(texts, origPos):
                        text._x = v
                        ax.plot([0., text._x], [opos, text._y], color=text._color, lw=0.5, clip_on=False)

                # Clusters outline boxes
                if not plottingMajorMetricOfSelected:
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

                    ltext = '' # 'Eucl. dist. of gene expr.\n %s dist.' % self.majorMetric
                    clb = fig.colorbar(im, ax=ax, fraction=0.4, label=ltext)
                    clb.ax.tick_params(labelsize=fontsize)

                return

            def addBar(fig, dataArgs, panel, coords, halfWindowSize = halfWindowSize, noAverage = False):

                nonlocal panelsData, panelsDataNames, externalPanelsData, togglePublicationFigure

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
                    if togglePublicationFigure:
                        ylabel='Expression\nfraction'
                    else:
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
                    if togglePublicationFigure:
                        ylabel='Known angiogenesis\nreceptors'
                    else:
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
                    if togglePublicationFigure:
                        ylabel='Evolutionary\nconservation'
                    else:
                        ylabel='Top50\noverlap'
                    try:
                        data = externalPanelsData['conservedGenes']
                        data = data.loc[~data.index.duplicated(keep='first')].reindex(allGenes).values
                    except:
                        data = np.zeros(len(allGenes))

                elif panel == 'binomial':
                    if togglePublicationFigure:
                        ylabel='Network\nenrichment'
                    else:
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

                    if ylabel == 'combo3avgs':
                        if togglePublicationFigure:
                            ylabel = 'Combination\nof measures' # 'Combination of 3'
                    elif ylabel == 'combo4avgs':
                        if togglePublicationFigure:
                            ylabel = 'Combination of 4'

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

                else:
                    ylabel = panel
                    try:
                        data = externalPanelsData[panel]
                        data = data.loc[~data.index.duplicated(keep='first')].reindex(allGenes).values
                    except:
                        data = np.zeros(len(allGenes))

                ax.bar(range(len(clusters)), data, width=ax.get_xlim()[1]/len(clusters), color=tickLabelsColors)

                data_avg = movingAverageCentered(np.nan_to_num(data), halfWindowSize)

                if not noAverage:
                    ax.plot(range(len(clusters)), data_avg, linewidth=1.0, color='coral', alpha=1.0)

                    if not togglePublicationFigure:
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

                ax.tick_params(axis='y', labelsize=4, width=0.75, length=3)
                ax.yaxis.tick_right()

                if showBoundaries:
                    ylim = ax.get_ylim()
                    for i in range(1, len(np.unique(clusters))):
                        ax.plot([clusterBoundaries[i] - 0.5]*2, [ylim[0], ylim[1]], '--', lw=0.5, color='k', clip_on=False)

                ax.text(-0.01, 0.5, ylabel, fontsize=6.5, rotation=0, va='center', ha='right', transform=ax.transAxes)
                
                if togglePublicationFigure:
                    if panel == 'combo3avgs':
                        peaks = getGenesOfPeak(pd.Series(data/data.max()))
                        bw = ax.get_xlim()[1]/len(clusters)
                        ax.axvspan(min(peaks) - 0.5*bw, max(peaks) + 0.5*bw, facecolor='crimson', linewidth=0., zorder=-np.inf)

                return ylabel, data, data_avg


            mmin, mmax = np.nanmin(np.nanmin(df.values)), np.nanmax(np.nanmax(df.values))

            if self.majorMetric in ['correlation', 'spearman', 'cosine']:
                missingFillValue = 1.0
            elif self.majorMetric == 'euclidean':
                missingFillValue = mmax
            else:
                missingFillValue = mmax

            if printStages:
                print('\tCalculating %s metric of %s . . .' % (metric, self.majorMetric), end='\t', flush=True)

            df = df.fillna(missingFillValue)

            if True:
                M = pdist(df.values.T, metric=metric)
            else:
                M = df.loc[df.columns].values
                M = (M + M.T)/2.
                np.fill_diagonal(M, 0.)
                M = squareform(M)

            if printStages:
                print('Done', flush=True)

            nPanels = len(self.panels)
            dendroHeight = 0.10
            panelHeight = 0.030 # 0.022
            delta = 0.01

            vsum = nPanels*(panelHeight + delta) + dendroHeight + 0.8

            if toggleAdjustFigureHeight:
                figureSize = (figureSize[0],  figureSize[0] * (9./8.) * vsum)

            factor = 1 / vsum

            topBorder = 1. - 0.03*factor
            bottomBorder = 0.05*factor

            dendroHeight *= factor
            panelHeight *= factor
            delta *= factor

            heatmapHeight = (topBorder - bottomBorder - dendroHeight) - nPanels * (panelHeight + delta) - 0.05*factor

            if not toggleIncludeHeatmap:
                rescaleFactor = 1. - heatmapHeight
                figureSize = figureSize[0], figureSize[1] * rescaleFactor
                heatmapHeight = 0.
                topBorder = 1. - (1. - topBorder) / rescaleFactor
                dendroHeight /= rescaleFactor
                panelHeight /= rescaleFactor
                delta /= rescaleFactor
                bottomBorder /= rescaleFactor

            fig = plt.figure(figsize=figureSize)

            dataArgs = addDendro(fig, df.columns, M, [0.1, topBorder-dendroHeight, 0.75, dendroHeight], metric=metric, linkageMethod=linkageMethod)

            if nPanels > 0:
                if printStages:
                    print('\tPlotting bar panels . . .', end='\t', flush=True)
                panelsData = dict()
                panelsDataNames = dict()
                for ipanel, panel in enumerate(self.panels):
                    coords = [0.1, 0.015*factor + bottomBorder + heatmapHeight + (len(self.panels) - ipanel - 1)*(panelHeight + delta), 0.75, panelHeight]

                    if panel[-1] == '~':
                        fig.add_artist(lines.Line2D([0, 0.875], [coords[1] - 0.5 * delta]*2, linewidth=0.75))
                        panel = panel[:-1]

                    panelName, data, data_avg = addBar(fig, dataArgs, panel, coords)

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
                    dataArgs.update({'majorMetricOfSelected': 1. - df.loc[dataArgs['allGenes']][dataArgs['allGenes']].values})
                    addHeatmap(fig, dataArgs, [0.1, bottomBorder, 0.75, heatmapHeight])
                    if printStages:
                        print('Done', flush=True)

                number = np.loadtxt(os.path.join(saveDir, 'size.txt'), dtype=int)
                np.savetxt(os.path.join(saveDir, 'figureInfo.txt'), np.array([['Cells', number[1]], ['Genes', number[0]], ['SelGenes', df.shape[1]]]), fmt='%s')

                if not togglePublicationFigure:
                    fig.suptitle('Ordering by single cell co-expression\nData: %s\n(%s cells, %s x %s genes)' % (suffix, number[1], number[0], df.shape[1]), fontsize=8, backgroundcolor='white')

                if printStages:
                    print('\tSaving image . . .', end='\t', flush=True)
                fig.savefig(os.path.join(saveDir, '%s dendrogram-heatmap-%s.png' % (suffix, self.majorMetric)), dpi=dpi)
                if printStages:
                    print('Done', flush=True)

            plt.close(fig)

            return dataArgs

        def exportFigureData(dataArgs, saveXLSX = True, saveHDF = True):

            '''Function to assist in exporting figure data to excel and/or hdf file 

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

    def reanalyzeMain(self, case='All', **kwargs):

        '''Reanalyze case 

        Parameters: 
            Any parameters that function 'analyzeCase' can accept

        Returns:
            None 

        Usage:
            an = Analyze()

            an.reanalyzeMain()
        '''

        np.random.seed(42)

        try:
            print('Re-analyzing %s-data case' % case, flush=True)

            kwargs.setdefault('toggleAdjustText', True)
            kwargs.setdefault('dpi', 600)
            kwargs.setdefault('suffix', case)
            kwargs.setdefault('saveDir', os.path.join(self.bootstrapDir, '%s/' % case))
            kwargs.setdefault('printStages', True)
            kwargs.setdefault('toggleCalculateMajorMetric', False)
            kwargs.setdefault('toggleExportFigureData', True)
            kwargs.setdefault('toggleCalculateMeasures', True)

            kwargs.setdefault('externalPanelsData', externalPanelsData)
            kwargs['externalPanelsData'].update({'conservedGenes': pd.read_excel(os.path.join(self.bootstrapDir, case, 'comparison.xlsx'), index_col=1, header=0)['Inter-measures.T50_common_count']})
            
            self.analyzeCase(None, **kwargs)

            shutil.copyfile(os.path.join(self.bootstrapDir, case, '%s dendrogram-heatmap-%s.png' % (case, self.majorMetric)), os.path.join(self.workingDir, 'results %s %s %s %s.png' % (case, self.majorMetric, self.dendrogramMetric, self.dendrogramLinkageMethod)))
            
            shutil.copyfile(os.path.join(self.bootstrapDir, case, 'dendrogram-heatmap-%s-data.xlsx' % (self.majorMetric)), os.path.join(self.workingDir, 'results %s %s %s %s.xlsx' % (case, self.majorMetric, self.dendrogramMetric, self.dendrogramLinkageMethod)))

        except Exception as exception:
            print(exception)

        return

    def analyzeAllPeaksOfCombinationVariant(self, variant, nG = 8, nE = 30, dcutoff = 0.5, fcutoff = 0.5, width = 50):

        '''Find all peaks and their frequency from the bootstrap experiments

        Parameters: 
            variant: str
                Name of combination variant (e.g. 'Avg combo4avgs', 'Avg combo3avgs')

            nG: int, Default 8
                Number of clusters of genes

            nE: int, Default 30
                Number of clusters of bootstrap experiments

            fcutoff: float, Default 0.5
                Lower peak height cutoff

            width: int, Default 50
                Width of peak 

        Returns:
            None 

        Usage:
            an = Analyze()
            
            an.analyzeAllPeaksOfCombinationVariant('Avg combo4avgs', nG=8, nE=30, fcutoff=0.5, width=50)
        '''

        def getPeaksLists(df_temp):

            '''Get list of peaks and genes

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
                peaksLists.update({(experiment, i, se.iloc[peak]): getGenesOfPeak(se, peak=peak, maxDistance=int(width/2)) for i, peak in enumerate(getPeaks(se, threshold = 0.05, prominence = 0.05, distance=width))})
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
        print(df_m.shape)
        df_m = df_m.loc[df_m.sum(axis=1) > (0.05 * df_m.shape[1])]
        print(df_m.shape)

        se_heights = pd.Series(index=['E' + df_m.columns.get_level_values(0).str.split('Experiment ', expand=True).get_level_values(-1) + '.' + df_m.columns.get_level_values(1).astype(str)], data=df_m.columns.get_level_values(-1))
        se_heights.index = se_heights.index.get_level_values(0)

        # Example plot of determining peaks
        if False:
            e = 'All' # 'Experiment 1' 'All'
            
            if False:
                self.reanalyzeMain(case=e, togglePublicationFigure=True, toggleIncludeHeatmap=False, markersLabelsRepelForce=1.5)
            
            if e == 'All':
                se = pd.read_hdf(os.path.join(self.bootstrapDir, 'All', 'dendrogram-heatmap-%s-data.h5' % self.majorMetric), key='df_C')
                #se = se[variant]
                se = se['Avg Combination']
                se.index.name = 'gene'
            else:
                se = df.xs('species', level='species', axis=0).copy().xs(e, level='experiment')

            fig, ax = plt.subplots(figsize=(8,2))
            ax.plot(range(len(se)), se.values, color='coral', linewidth=1.5, zorder=np.inf)

            if e == 'All':
                listsDict, allgenes = getPeaksLists(pd.concat([se.to_frame().copy()], keys=['All'], names=['experiment'], axis=0, sort=False)['Avg Combination'])
                #listsDict, allgenes = getPeaksLists(pd.concat([se.to_frame().copy()], keys=['All'], names=['experiment'], axis=0, sort=False)[variant])
                se_peakAssignments = pd.Series(listsDict).apply(pd.Series)
                df_m_a = pd.DataFrame(index=allgenes, data=0, columns=se_peakAssignments.index)
                for peak in df_m_a.columns:
                    df_m_a.loc[se_peakAssignments.loc[peak].dropna().values, peak] = 1
                se_peaks = df_m_a.xs(e, level=0, axis=1)
            else:
                se_peaks = df_m.xs(e, level=0, axis=1)

            for peak in se_peaks[:]:
                se_peak = se_peaks[peak]
                se_peak = se_peak[se_peak == 1]
                x = np.where(np.isin(se.index.values, se_peak.index.values))[0]
                y = se.iloc[x].values

                if peak[1] == 1.:
                    print(len(se_peak))
                    print(cleanListString(se_peak.index.values.tolist()))

                c = 'grey' if peak[1] != 1. else 'blue'
                ax.fill_between(x, y, facecolor=c, edgecolor=c, alpha=0.5)

                pg = se.index.values[x][np.argmax(y)]
                ax.text(x[np.argmax(y)], -0.001, pg, fontsize=8, rotation=90, va='top', ha='center')
                ax.plot([x[np.argmax(y)] - 1, x[np.argmax(y)] - 1], [0, max(y)], color='k', linewidth=1., alpha=1.)

            ax.set_xticks([])
            ax.set_xticklabels([])
            yticks = (0, 0.03)
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticks)
            ax.tick_params(axis='y', labelsize=8, width=0.75, length=3)
            ax.text(-0.05, 0.5, 'Combination\nof measures', fontsize=10, rotation=90, va='center', ha='right', transform=ax.transAxes)
            ax.set_xlim([0, len(se)])
            ax.set_ylim([0, 0.03])

            fig.tight_layout()
            fig.savefig(self.workingDir + 'one example peaks %s.png' % variant, dpi=600)

            plt.close(fig)

            exit()

        Z1 = hierarchy.linkage(df_m, 'ward')
        Z2 = hierarchy.linkage(df_m.T, 'ward')

        O1 = hierarchy.dendrogram(Z1, no_plot=True, get_leaves=True)['leaves'][::-1]
        O2 = hierarchy.dendrogram(Z2, no_plot=True, get_leaves=True)['leaves']

        df_m = df_m.iloc[O1]
        df_m = df_m.T.iloc[O2].T

        clusters1 = hierarchy.fcluster(Z1, t=nG, criterion='maxclust')[O1] - 1
        clusters2 = hierarchy.fcluster(Z2, t=nE, criterion='maxclust')[O2] - 1

        u1 = np.unique(clusters1, return_counts=True)[1]
        u2 = np.unique(clusters2, return_counts=True)[1]

        df_m.index = pd.MultiIndex.from_arrays([df_m.index, clusters1], names=['gene', 'cluster'])
        df_m.columns = pd.MultiIndex.from_arrays(['E' + df_m.columns.get_level_values(0).str.split('Experiment ', expand=True).get_level_values(-1) + '.' + df_m.columns.get_level_values(1).astype(str), clusters2], names=['peak', 'cluster'])

        we = pd.Series(index=df_m.columns, data=se_heights[df_m.columns.get_level_values(0)].values).droplevel(1) # 570 peaks
        ndict = dict()
        resf = dict()
        resh = dict()
        for ci in df_m.index.levels[-1]:
            for cj in df_m.columns.levels[-1]:
                df_temp = df_m.xs(ci, level='cluster', axis=0).xs(cj, level='cluster', axis=1)
                if df_temp.values.mean() >= dcutoff:
                    resf[(ci, cj)] = df_temp.shape[1]
                    resh[(ci, cj)] = we[df_temp.columns].mean()

                    ndict[ci] = cleanListString(sorted(df_temp.index.values.tolist()))

        sef = pd.Series(resf)
        sed = pd.Series(index=sef.index, data=sef.groupby(level=0).sum().reindex(sef.droplevel(1).index).values)
        dff = pd.concat([pd.Series(resh), sef / sed], axis=1, keys=['h', 'w'], sort=False)
        seh = dff.groupby(level=0).agg(lambda s: np.average(s, weights=dff.loc[s.index, 'w'])).h

        #pd.Series(resf).to_excel(self.workingDir + 'resf %s.xlsx' % variant)
        #pd.Series(resh).to_excel(self.workingDir + 'resh %s.xlsx' % variant)

        se = pd.Series(resf).groupby(level=0).sum()
        se.name = 'frequency'
        df = se.to_frame()
        df['height'] = seh
        df['genes'] = pd.Series(se.index).replace(ndict).values
        df['length'] = df['genes'].str.split(', ').map(len)
        df['frequency'] = np.zeros(len(df['frequency']))
        df = df[['frequency', 'height', 'length', 'genes']]

        locnbootstap = np.unique(se_peakAssignments.index.get_level_values(0)).shape[0]

        for i, group in enumerate(df['genes'].values):
            group = group.split(', ')
            fractions = df_m.droplevel('cluster').loc[group].sum(axis=0) / len(group)
            uv = np.unique(fractions[fractions >= fcutoff].index.get_level_values(0).str.split('.', expand=True).get_level_values(0).values)
            df.loc[df.index[i], 'frequency'] = len(uv) / locnbootstap

        df = df.sort_values(by='frequency', ascending=False)
        print(df)

        df.columns = ['Fequency of group', 'Average peak height', 'Number of genes', 'Genes']
        df.index = np.array(range(len(df.index))) + 1
        df.index.name = 'Group'
        df.to_excel(self.workingDir + 'All peaks %s.xlsx' % variant)

        # Enrichment analysis
        if False:
            from pyiomica.enrichmentAnalyses import KEGGAnalysis, GOAnalysis, ReactomeAnalysis, ExportEnrichmentReport, ExportReactomeEnrichmentReport

            dataForAnalysis = df['Genes'].str.split(', ').to_dict()

            resKEGG = KEGGAnalysis(dataForAnalysis)
            ExportEnrichmentReport(resKEGG, AppendString='KEGG11', OutputDirectory=self.workingDir)

            resGO = GOAnalysis(dataForAnalysis)
            for key in resGO.keys():
                resGO[key] = {k:v for k,v in resGO[key].items() if v[2][1]=='biological_process'}
            ExportEnrichmentReport(resGO, AppendString='GO11', OutputDirectory=self.workingDir)
            
            resReactome = ReactomeAnalysis(dataForAnalysis)
            resReactome = {k:v.loc[v['Entities FDR'] < 0.05] for k,v in resReactome.items()}
            ExportReactomeEnrichmentReport(resReactome, AppendString='Reactome11', OutputDirectory=self.workingDir)

        # Plot all-peaks heatmap
        if True:
            fig = plt.figure(figsize=(10, 10))

            groupsColors = True

            # Genes dendrogram
            if True:
                n_clusters = nG
                ax = fig.add_axes([0.1, 0.25, 0.15, 0.5], frame_on=False)
                origLineWidth = matplotlib.rcParams['lines.linewidth']
                matplotlib.rcParams['lines.linewidth'] = 0.5
                if groupsColors:
                    D = hierarchy.dendrogram(Z1, ax=ax, color_threshold=Z1[-n_clusters + 1][2], above_threshold_color='k', orientation='left')
                else:
                    D = hierarchy.dendrogram(Z1, ax=ax, color_threshold=0, above_threshold_color='k', orientation='left')
                matplotlib.rcParams['lines.linewidth'] = origLineWidth
                ax.set_xticklabels([])
                ax.set_xticks([])
                ax.set_yticklabels([])
                ax.set_yticks([])
                ax.set_title('Receptors')

            # Peaks dendrogram
            if True:
                n_clusters = nE
                ax = fig.add_axes([0.25, 0.75, 0.5, 0.15], frame_on=False)
                origLineWidth = matplotlib.rcParams['lines.linewidth']
                matplotlib.rcParams['lines.linewidth'] = 0.5
                if groupsColors:
                    D = hierarchy.dendrogram(Z2, ax=ax, color_threshold=Z2[-n_clusters + 1][2], above_threshold_color='k', orientation='top')
                else:
                    D = hierarchy.dendrogram(Z2, ax=ax, color_threshold=0, above_threshold_color='k', orientation='top')
                matplotlib.rcParams['lines.linewidth'] = origLineWidth
                ax.set_xticklabels([])
                ax.set_xticks([])
                ax.set_yticklabels([])
                ax.set_yticks([])
                ax.set_title('Peaks')

            # Heatmap
            if True:
                heatdata = df_m.values * we.values[None, :]
                ax = fig.add_axes([0.25, 0.25, 0.5, 0.5], frame_on=True)
                #cmap = matplotlib.colors.LinearSegmentedColormap.from_list('WB', [(1, 1, 1), (0, 0, 0.5)], N=2)
                cmap = copy.copy(plt.cm.jet_r)
                cmap.set_bad('white')
                im = ax.imshow(np.ma.array(heatdata, mask=(heatdata==0.)), cmap=cmap, aspect='auto', interpolation='None', extent=(-0.5, df_m.shape[0] - 0.5, df_m.shape[1] - 0.5, -0.5), vmin=0., vmax=1.)
                ax.set_xticklabels([])
                ax.set_xticks([])
                ax.set_yticklabels([])
                ax.set_yticks([])

            # Export heatmap data
            if True:
                print('Exporting all peaks heatmap', flush=True)
                (df_m * we.values[None, :]).to_excel(self.workingDir + 'heatmap all peaks %s.xlsx' % variant)
            
            # Box annotations
            if False:
                xl, yl = ax.get_xlim()[1], ax.get_ylim()[0]
                sea = pd.Series(index=pd.MultiIndex.from_tuples(df_m.index.values), data=df_m.index.values).apply(pd.Series)
                gsea = sea.groupby(level=1).count()[0]
                dexp = pd.Series(index=df_m.columns.get_level_values(1).values, data=df_m.columns.get_level_values(1).values)
                gexp = dexp.groupby(level=0).count()

                def ann(xy2, clusterG, clusterE, n = 10):

                    xy1 = xl * np.mean(np.where(dexp == clusterE)[0]) / len(dexp), yl * np.mean(np.where(sea[1] == clusterG)[0]) / len(sea)
                    xy2 = xy2[0] * xl, xy2[1] * yl

                    genes = sea.xs(key=clusterG, level=1)[0].values
                    sgenes = '\n'.join([cleanListString(g) for g in np.array_split(genes, int(len(genes)/n))])

                    ax.annotate(sgenes, xy=xy1, xytext=xy2, fontsize=5, ha='left', va='center', arrowprops=dict(arrowstyle="-", lw=0.5), bbox=dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.25), zorder=np.inf)

                    return
                
                ann((0.76, 1.14), 8, 6, n=5)
                ann((0.76, 1.05), 7, 7, n=5)

                ann((0.39, 1.09), 1, 2, n=5)

                ann((-0.3, 1.125), 0, 0, n=10)
                ann((-0.3, 1.05), 6, 1, n=10)

            # Colorbar
            if True:
                axColor = fig.add_axes([0.72, 0.25, 0.05, 0.3], frame_on=False)
                axColor.set_xticks([])
                axColor.set_xticklabels([])
                axColor.set_yticks([])
                axColor.set_yticklabels([])
                clb = fig.colorbar(im, ax=axColor)
                clb.set_label(label='Peak height', fontsize=6)
                #clb.set_ticks([0.25, 0.75])
                #clb.set_ticklabels(['No', 'Yes'])
                clb.ax.tick_params(labelsize=6)
            
            fig.savefig(self.workingDir + 'all peaks %s.png' % variant, dpi=600)
            plt.close(fig)

        return
        
    def analyzePerGeneCombinationVariant(self, variant, hcutoff = 0.2, fcutoff = 0.3, width = 50):

        '''Find all peaks and their frequency from the bootstrap experiments

        Parameters: 
            variant: str
                Name of combination variant (e.g. 'Avg combo4avgs', 'Avg combo3avgs')

            nG: int, Default 8
                Number of clusters of genes

            nE: int, Default 30
                Number of clusters of bootstrap experiments

            fcutoff: float, Default 0.5
                Lower peak height cutoff

            width: int, Default 50
                Width of peak 

        Returns:
            None 

        Usage:
            an = Analyze()
            
            an.analyzeAllPeaksOfCombinationVariant('Avg combo4avgs', nG=8, nE=30, fcutoff=0.5, width=50)
        '''

        print('Variant:', variant)
        df = pd.read_hdf(self.dendroDataName, key='df').fillna(0).set_index('gene', append=True).droplevel('order')[variant]

        df_temp = df.xs('species', level='species', axis=0).copy()
        experiments = np.unique(df_temp.index.get_level_values('experiment'))
        genes = np.unique(df_temp.index.get_level_values('gene'))

        dfs = []
        for gene in genes:
            peaksLists = {}
            for experiment in experiments:
                se = df_temp.xs(experiment, level='experiment', axis=0)
                se = se/se.max()

                if gene in se.index:
                    i0 = np.where(se.index==gene)[0][0]
                    i1, i2 = max(0, i0-int(width/2)), min(i0+int(width/2), len(se))

                    # e.g. w=50, then take 25 to the left of the gene, the gene itself, 24 to the right of the gene
                    se_temp = se.iloc[i1: i2]
                    se_temp = se_temp[se_temp >= hcutoff]
                    peaksLists.update({experiment: se_temp.index.values})

            df_peaks = pd.Series(peaksLists).apply(pd.Series)
            se = df_peaks.stack().dropna().value_counts().sort_values(ascending=False)/len(experiments)
            se = se.reindex(genes).fillna(0.)
            se.name = gene

            if se[gene] >= fcutoff:
                print(gene, end=' ', flush=True)
                dfs.append(se)

        dfs = pd.concat(dfs, axis=1, sort=False)
        print(dfs)
        dfs.to_excel(self.workingDir + 'near frequency %s %s.xlsx' % (variant, hcutoff))

        return

    def bootstrapMaxpeakPlot(self, variant):

        '''Bootstrap max-peak plot
        '''

        df = pd.read_excel(self.workingDir + '%s bootstrap_in-peak_genes_SD.xlsx' % variant, sheet_name='Bootstrap lists', index_col=None, header=0)
        df_temp = pd.DataFrame(index=np.unique(df.stack().dropna().values), columns=df.columns).fillna(0.)
        for col in df.columns:
            df_temp.loc[df[col].dropna().values, col] = 1.

        df_temp = df_temp.iloc[scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(df_temp, 'ward'), no_plot=True, get_leaves=True)['leaves']]
        df_temp = df_temp.T.iloc[scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(df_temp.T, 'ward'), no_plot=True, get_leaves=True)['leaves']].T

        df_temp = df_temp.loc[df_temp.sum(axis=1) > (0.1 * df_temp.shape[1])]

        fig = plt.figure(figsize=(12, 10))

        ax = fig.add_axes([0.25, 0.2, 0.65, 0.7])
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list('WB', [(1, 1, 1), (0, 0, 0.5)], N=2)
        im = ax.imshow(df_temp.values, cmap=cmap, interpolation='None', aspect='auto')

        ax.set_xticks(range(df_temp.shape[1]))
        ax.set_xticklabels(df_temp.columns, rotation=90, fontsize=5)
        ax.set_xlim([-0.5, df_temp.shape[1] - 0.5])

        ax.set_yticks(range(df_temp.shape[0]))
        ax.set_yticklabels(df_temp.index, rotation=0, fontsize=5)
        ax.set_ylim([-0.5, df_temp.shape[0] - 0.5])

        if True:
            axColor = fig.add_axes([0.825, 0.8, 0.1, 0.1], frame_on=False)
            axColor.set_xticks([])
            axColor.set_xticklabels([])
            axColor.set_yticks([])
            axColor.set_yticklabels([])
            clb = fig.colorbar(im, ax=axColor)
            clb.set_ticks([0.25, 0.75])
            clb.set_ticklabels(['No', 'Yes'])
            clb.ax.tick_params(labelsize=6)

        fig.savefig(self.workingDir + '%s bootstrap.png' % variant, dpi=600)

        return

    def generateAnalysisReport(self):

        '''Generate analysis report.
        '''
        
        saveDir = os.path.join(self.workingDir, 'results majorMetric=%s, dendrogramMetric=%s, linkageMethod=%s' % (self.majorMetric, self.dendrogramMetric, self.dendrogramLinkageMethod))

        if not os.path.exists(saveDir):
            os.makedirs(saveDir)

        for file in os.listdir(self.workingDir):
            if file != 'data.h5' and not os.path.isdir(os.path.join(self.workingDir, file)):
                try:
                    shutil.copyfile(os.path.join(self.workingDir, file), os.path.join(saveDir, file))
                except Exception as exception:
                    print('Cannot copy file:', file, '\n', exception)

        return 

def process(df1main, df1other, df2main, df2other, dir1, dir2, genesOfInterest = None, knownRegulators = None, nCPUs = 4,
            panels = ['fraction', 'binomial', 'top50', 'markers', 'combo3avgs', 'combo4avgs'], parallelBootstrap = False,
            exprCutoff1 = 0.05, exprCutoff2 = 0.05, perEachOtherCase = True, doScramble = False, part1 = True, part2 = True,
            part3 = True, **kwargs):

    '''Main workflow programmed in two scenaria depending on parameter "perEachOtherCase".

    Parameters:
        df1main: pandas.DataFrame
            Expression data of main group of cells of the first species

        df1other: pandas.DataFrame
            Expression data of other cells of the first species

        df2main: pandas.DataFrame
            Expression data of main group of cells of the second species

        df2other: pandas.DataFrame
            Expression data of other cells of the second species

        dir1: str
            Path to the first species working directory

        dir2: str
            Path to the second species working directory

        genesOfInterest: list, Default None
            Particular genes to analyze, e.g. receptors

        knownRegulators: list, Default None
            Known marker genes 

        nCPUs: int, Default 1
            Number of CPUs to use for multiprocessing, recommended 10-20

        panels: list, Default None
            Particular measurements to include in the analysis

        parallelBootstrap: boolean, Default False
            Whether to generate bootstrap experiments in parallel mode

        exprCutoff1: float, Default 0.05
            Per-batch expression cutoff for the first dataset

        exprCutoff2: float, Default 0.05
            Per-batch expression cutoff for the second dataset

        perEachOtherCase:  boolean, Default True
            Scenario of comparison

        Any other parameters that class "Analysis" can take

    Returns:
        Analysis
            First class Analysis instance

        Analysis
            Second class Analysis instance
    '''

    if genesOfInterest is None:
        genesOfInterest = receptorsListHugo_2555
    
    if knownRegulators is None:
        knownRegulators = gEC22

    an1 = Analysis(workingDir=dir1, otherCaseDir=dir2, genesOfInterest=genesOfInterest,
                        knownRegulators=knownRegulators, panels=panels, nCPUs=nCPUs, perEachOtherCase=perEachOtherCase, **kwargs)

    if perEachOtherCase:
        an2 = Analysis(workingDir=dir2, otherCaseDir=dir1, genesOfInterest=genesOfInterest,
                            knownRegulators=knownRegulators, panels=panels, nCPUs=nCPUs, perEachOtherCase=perEachOtherCase, **kwargs)

    if part1:
        if not df1main is None:
            an1.prepareDEG(df1main, df1other)

        an1.preparePerBatchCase(exprCutoff=exprCutoff1)
        an1.prepareBootstrapExperiments(parallel=parallelBootstrap)

        if perEachOtherCase:
            if not df2main is None:
                an2.prepareDEG(df2main, df2other)

            an2.preparePerBatchCase(exprCutoff=exprCutoff2)
            an2.prepareBootstrapExperiments(parallel=parallelBootstrap)

    if part2:
        an1.analyzeBootstrapExperiments()

        if perEachOtherCase:
            an2.analyzeBootstrapExperiments()
            an1.analyzeBootstrapExperiments()
            
        an1.reanalyzeMain()

        if an1.nBootstrap > 0:
            an1.analyzeCombinationVariant('Avg combo3avgs')
            an1.analyzeCombinationVariant('Avg combo4avgs')
            an1.bootstrapMaxpeakPlot('Avg combo3avgs')
            an1.bootstrapMaxpeakPlot('Avg combo4avgs')
            an1.analyzeAllPeaksOfCombinationVariant('Avg combo3avgs', nG=15, nE=15, fcutoff=0.5, width=50)
            an1.analyzeAllPeaksOfCombinationVariant('Avg combo4avgs', nG=15, nE=15, fcutoff=0.5, width=50)

        if doScramble:
            an1.scramble(['Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo3/', M=20)
            an1.scramble(['Markers', 'Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo4/', M=20)

        if perEachOtherCase:
            an2.reanalyzeMain()

            if an2.nBootstrap > 0:
                an2.analyzeCombinationVariant('Avg combo3avgs')
                an2.analyzeCombinationVariant('Avg combo4avgs')
                an2.bootstrapMaxpeakPlot('Avg combo3avgs')
                an2.bootstrapMaxpeakPlot('Avg combo4avgs')
                an2.analyzeAllPeaksOfCombinationVariant('Avg combo3avgs', nG=15, nE=15, fcutoff=0.5, width=50)
                an2.analyzeAllPeaksOfCombinationVariant('Avg combo4avgs', nG=15, nE=15, fcutoff=0.5, width=50)

            if doScramble:
                an2.scramble(['Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo3/', M=10)
                an2.scramble(['Markers', 'Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo4/', M=20)

    if part3:
        an1.generateAnalysisReport()

        if perEachOtherCase:
            an2.generateAnalysisReport()

    if perEachOtherCase:
        return an1, an2

    return an1
