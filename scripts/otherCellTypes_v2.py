from decneo.commonFunctions import *
#from decneo.supplementaryFunctions import getNonParametricPValue, getLogisticRegressionPValue
from decneo.analysisPipeline import Analysis, process
from decneo.genes import *


workingDir = '/mnt/home/paterno1/Ben/otherCellTypes/results/'
processedDataDir = '/mnt/research/piermarolab/Sergii/processedDataDir/'
bigDataPath = processedDataDir + 'selectedCellsInPanglaoDB_processed_4'
selCellsFile = processedDataDir + 'selectedCellsInPanglaoDB_4.h5'

if __name__ == '__main__':

    # Count SRA/SRS to find cell types with plenty (>20) SRS per both, human and mouse
    # Add subtypes of the selected cell types too, name conversion is in "dictCelltypes"
    if False:
        df = pd.read_excel('dev/PanglaoDBdata/df_cell_type_annotations.xlsx', index_col=[0, 1], header=0)[['Cell type annotation', 'Species']].reset_index()[['Species', 'Cell type annotation', 'SRA accession']].drop_duplicates()
        df = df.set_index(['Species', 'Cell type annotation'])['SRA accession']
        print(df)

        dfg = df.groupby(level=['Species', 'Cell type annotation'], axis=0).count().unstack('Species').sort_values('Mus musculus', ascending=False).reset_index()
        print(dfg)
        dfg.to_excel('dfgSRA.xlsx')

    #uCelltypes = ['Pericyte', 'SMC']
    #dictCelltypes = {'Pericytes': 'Pericyte', 'Smooth muscle cells': 'SMC'}

    uCelltypes = ['Macrophage']
    dictCelltypes = {'Macrophages': 'Macrophage'}

    # Load all selected cells IDs into one dataframe
    if False:
        df_cell_type_annotations = getPanglaoDBAnnotationsSummaryDf(MetadataDirName)[['Cell type annotation', 'Species']]
        df_cell_type_annotations = df_cell_type_annotations.loc[np.isin(df_cell_type_annotations['Cell type annotation'].values, list(dictCelltypes.keys()))]
        df_cell_type_annotations['Cell type annotation'] = df_cell_type_annotations['Cell type annotation'].replace(dictCelltypes)
        df_cell_type_annotations = df_cell_type_annotations.set_index(['Cell type annotation', 'Species'], append=True)
        print(df_cell_type_annotations.shape)
    
        dfs = []
        for SRA, SRS, cluster, celltype, species in df_cell_type_annotations.index[:]:
            df_cell_to_cluster = pd.read_csv(os.path.join(MetadataDirName, 'PanglaoDB', 'data', 'sample_clusters', '%s%s.seurat_clusters.txt' % (SRA, '_' + SRS if SRS!='notused' else '')), delimiter=' ', index_col=0, header=None)[1]
            selectedCells = df_cell_to_cluster[df_cell_to_cluster == cluster].index.values
            dfs.append(pd.DataFrame(index=pd.MultiIndex.from_tuples([(species, '%s%s.sparse.RData.h5' % (SRA, '_' + SRS if SRS!='notused' else ''))], names=['species', 'file']), columns=np.unique(selectedCells), data=celltype).stack())

        dfs = pd.concat(dfs, axis=0, sort=False)
        dfs.to_hdf(selCellsFile, key='df', mode='a', complevel=4, complib='zlib')
        print(dfs)

    # Counts summary
    if False:
        df = pd.read_hdf(selCellsFile, key='df').to_frame()
        df.index.names = ['species', 'file', 'cell']
        df.columns = ['celltype']

        df = df.reset_index().set_index(['species', 'celltype', 'file'])['cell']
        df = df.groupby(['species', 'celltype', 'file']).count().unstack('celltype')
        print(df)

        df.reset_index().to_excel(selCellsFile + 'counts.xlsx', index=False)

    # Prepare expression and DEG ranking
    if False:
        def parseFiles(filesP):

            df_cells_id = pd.read_hdf(selCellsFile, key='df')
            for species in ['Homo sapiens', 'Mus musculus'][:]:
                df_cells_id_species = df_cells_id.xs(species, level='species')
                files = np.unique(df_cells_id_species.index.get_level_values('file'))
                print('Species:', species, '\tFiles:', len(files), '\tCells:', len(df_cells_id_species), flush=True)

                for ifile, file in enumerate(files):
                    if file not in filesP:
                        continue

                    saveFileNamePath = bigDataPath + '/%s.h5' % file

                    if os.path.isfile(saveFileNamePath +'.txt'):
                        continue

                    try:
                        print('\nProcessing file:', ifile, file, 'of', len(files), species, flush=True)

                        # Load expression data
                        df_expr = pd.read_hdf(os.path.join(RDataDirName, file), key='df')

                        # Make sure all cell identifiers are found in expression file
                        df_temp_cells = df_cells_id_species.xs(file, level='file')
                        df_temp_cells = df_temp_cells.loc[np.isin(df_temp_cells.index.values, df_expr.columns)]

                        # Use selected cells
                        w = df_expr.columns.difference(df_temp_cells.index)
                        df_expr = df_expr[pd.Index(np.random.choice(w, min(1000, len(w)), replace=False)).append(df_temp_cells.index)]

                        # Convert any mouse genes to human
                        df_expr.index = pd.Series(df_expr.index.values).replace(Mouse_to_Human_HUGO_conversion).values

                        # Drop any duplicates
                        df_expr = df_expr.loc[~df_expr.index.duplicated(keep='first')]
                        df_expr = df_expr.T.loc[~df_expr.T.index.duplicated(keep='first')].T
                        df_expr = df_expr.astype(float).fillna(0.)

                        # Drop zero cells
                        df_expr = df_expr[df_expr.columns[df_expr.sum(axis=0) > 0.]]

                        # Scale and log-transform
                        df_expr /= df_expr.sum(axis=0) * 0.0001
                        df_expr = np.log2(df_expr.replace(0., np.min(df_expr.values[df_expr.values > 0.])))
                        df_expr -= np.min(df_expr.values)

                        # Remove nearly-constant and constant genes
                        df_expr = df_expr[np.std(df_expr, axis=1) / np.mean(np.std(df_expr.values)) > 0.01]

                        for celltype in np.unique(df_temp_cells.values):
                            try:
                                celltypeCells = df_temp_cells[df_temp_cells == celltype].index
                                df_celltype = df_expr[celltypeCells]
                
                                # Rank and save differentially expressed genes
                                df_ttest = pd.DataFrame(index=df_celltype.index, columns=['statistic', 'pvalue'])
                                ttest = scipy.stats.ttest_ind(df_celltype.values, df_expr[df_expr.columns.difference(celltypeCells)].values, axis=1)
                                df_ttest['statistic'] = ttest[0]
                                df_ttest['pvalue'] = ttest[1]
                                df_ttest = df_ttest.sort_values('statistic', ascending=False).dropna()
                                df_ttest.to_hdf(saveFileNamePath, key='ttest/%s/%s/%s/' % (species, celltype, file), mode='a', complevel=4, complib='zlib')

                                # Save expression data
                                df_celltype = df_celltype[df_celltype.sum(axis=1) > 0]
                                df_celltype.to_hdf(saveFileNamePath, key='expr/%s/%s/%s/' % (species, celltype, file), mode='a', complevel=4, complib='zlib')

                                #print('\t', celltype, df_ttest.shape, df_celltype.shape, flush=True)

                            except Exception as exception:
                                print('Celltype ERROR:', exception, flush=True)

                        np.savetxt(saveFileNamePath +'.txt', ['Checked'], fmt='%s')

                    except Exception as exception:
                        print('ERROR:', exception, flush=True)

            return

        filesP = np.unique(pd.read_hdf(selCellsFile, key='df').index.get_level_values(1).values)
        print(filesP.shape)

        pool = multiprocessing.Pool(processes=36)
        pool.map(parseFiles, [[f] for f in filesP])
        pool.close()
        pool.join()

    # Aggregate files
    if False:
        keys = []
        allFiles = [file for file in os.listdir(bigDataPath) if file[-3:]=='.h5']
        for file in allFiles:
            keys.extend(KeysOfStore(bigDataPath + '/%s' % file))
        
        skeys = np.array([np.array(key.strip('/').split('/')) for key in keys])
        dfc = pd.DataFrame(skeys).set_index([0,1,2])
        dfc.index.names = ['x', 'species', 'celltype']
        dfc.columns = ['file']
        dfc = dfc.sort_index()
        print(dfc)

        for x in ['ttest', 'expr'][:]:
            for species in ['Homo sapiens', 'Mus musculus']:
                for celltype in uCelltypes:
                    try:
                        tempFiles = dfc.xs(x, level='x').xs(species, level='species').xs(celltype).values.flatten().tolist()

                        if x == 'ttest':
                            genes = []
                            batches = []
                            for file in tempFiles:
                                saveFileNamePath = bigDataPath + '/%s.h5' % file
                                df_temp = pd.read_hdf(saveFileNamePath, key='%s/%s/%s/%s/' % (x, species, celltype, file))
                                genes.append(df_temp.loc[df_temp['pvalue'] <= 10**-3]['statistic'].index.values)
                                batches.append(file.split('.sparse.RData.h5')[0])

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

                            print(df_ranks.head(5))
                            print('Saving ranks data:', species, celltype, df_ranks.shape, flush=True)
                            df_ranks.to_hdf(processedDataDir + 'PanglaoDB_ttest_ranks_per_batch_%s_%s.h5' % (species, celltype), key='df', mode='a', complevel=4, complib='zlib')

                        elif x == 'expr':
                            if os.path.isfile(processedDataDir + 'PanglaoDB_expresion_per_batch_%s_%s.h5' % (species, celltype)):
                                continue

                            dfs = []
                            for file in tempFiles:
                                batch = file.split('.sparse.RData.h5')[0]
                                saveFileNamePath = bigDataPath + '/%s.h5' % file
                                df_temp = pd.read_hdf(saveFileNamePath, key='%s/%s/%s/%s/' % (x, species, celltype, file))
                                df_temp = pd.concat([df_temp], keys=[batch], axis=1, sort=False)
                                df_temp.columns.names = ['batch', 'cell']
                                dfs.append(df_temp)

                            print('\nMerging expression data:', species, celltype, flush=True)
                            df_expr = pd.concat(dfs, axis=1, sort=False).fillna(0.)
                            
                            print('Saving expression data:', species, celltype, df_expr.shape, flush=True)
                            df_expr.to_hdf(processedDataDir + 'PanglaoDB_expresion_per_batch_%s_%s.h5' % (species, celltype), key='df', mode='a', complevel=4, complib='zlib')

                    except Exception as exception:
                        print('\nERROR:', exception, '\n')

    # Analyze per celltype
    if False:
        #mode = 'receptors' # 'ligands', 'receptors'

        try:
            mode = sys.argv[1]
        except Exception as exception:
            print('ARG ERROR:', exception)

        if mode=='ligands':
            knownRegulators = ligands_44
            genesOfInterest = ligands_1777
            workingDir = 'otherCellTypes_lig/'
        elif mode=='receptors':
            knownRegulators = gEC22
            genesOfInterest = receptorsListHugo_2555
            workingDir = 'otherCellTypes_rec/'
        else:
            raise NotImplementedError

        print('mode:', mode, '\t', 'knownRegulators:', len(knownRegulators), '\t', 'genesOfInterest:', len(genesOfInterest), flush=True)

        uCelltypes = ['Fibroblast', 'Endothelial']# Original uCelltypes = ['Macrophage', 'Pericyte', 'SMC', 'Fibroblast', 'Endothelial']
        print(uCelltypes)

        for celltype in uCelltypes:
            try:
                print(celltype)

                if True:
                    for species in ['Homo sapiens', 'Mus musculus']: 
                        an = Analysis(workingDir = workingDir + '%s/%s/' % (celltype, species), otherCaseDir = workingDir + '%s/%s/' % (celltype, 'Homo sapiens' if species=='Mus musculus' else 'Mus musculus'))

                        print('Saving DEG and expression data, %s, %s' % (species, celltype), flush=True)
                        df_expr = pd.read_hdf(processedDataDir + 'PanglaoDB_expresion_per_batch_%s_%s.h5' % (species, celltype), key='df')
                        df_expr.to_hdf(an.dataSaveName, key='df', mode='a', complevel=4, complib='zlib')
                        df_ranks = pd.read_hdf(processedDataDir + 'PanglaoDB_ttest_ranks_per_batch_%s_%s.h5' % (species, celltype), key='df')
                        df_ranks.to_hdf(an.dataSaveName, key='df_ranks', mode='a', complevel=4, complib='zlib')

                anHuman, anMouse = process(*(None, None), *(None, None), workingDir + '%s/Homo sapiens/' % celltype, workingDir + '%s/Mus musculus/' % celltype, nCPUs=20, parallelBootstrap=False, nBootstrap=100, PCNpath='/mnt/research/piermarolab/Sergii/results/', exprCutoff1=0.05, exprCutoff2=0.05, genesOfInterest=genesOfInterest, knownRegulators=knownRegulators, perEachOtherCase=True, part1=True, part2=True, part3=True)
    
            except Exception as exception:
                print('\nANALYSIS ERROR:', exception, '\n')

    # Re-plot main lig
    if True:
        for workingDir in ['otherCellTypes_rec/', 'otherCellTypes_lig/']:
            for celltype in ['Pericyte', 'SMC', 'Macrophage', 'Fibroblast', 'Endothelial']: # Original ['Pericyte', 'SMC', 'Macrophage', 'Fibroblast', 'Endothelial']:
                try:
                    if workingDir=='otherCellTypes_lig/':
                        knownRegulators = ligands_44
                        genesOfInterest = ligands_1777
                    elif workingDir=='otherCellTypes_rec/':
                        knownRegulators = gEC22
                        genesOfInterest = receptorsListHugo_2555

                    print(workingDir, celltype, len(knownRegulators), len(genesOfInterest))

                    anHuman, anMouse = process(*(None, None), *(None, None), workingDir + '%s/Homo sapiens/' % celltype, workingDir + '%s/Mus musculus/' % celltype, nCPUs=1, parallelBootstrap=False, nBootstrap=100, PCNpath='/mnt/research/piermarolab/Sergii/results/', exprCutoff1=0.05, exprCutoff2=0.05, panels=['fraction', 'binomial', 'top50', 'combo3avgs'], genesOfInterest=genesOfInterest, knownRegulators=knownRegulators, perEachOtherCase=True, part1=False, part2=False, part3=False)

                    anMouse.reanalyzeMain(togglePublicationFigure=True, includeClusterNumber=False, toggleIncludeHeatmap=True, toggleCalculateMeasures=True, toggleExportFigureData=True)
                    anHuman.reanalyzeMain(togglePublicationFigure=True, includeClusterNumber=False, toggleIncludeHeatmap=True, toggleCalculateMeasures=True, toggleExportFigureData=True)
                    try:
                        anHuman.analyzePerGeneCombinationVariant('Avg combo3avgs', hcutoff=0.2, fcutoff=0.0, width=50)
                    except Exception as exception:
                        print(exception)

                    try:
                        anMouse.analyzePerGeneCombinationVariant('Avg combo3avgs', hcutoff=0.2, fcutoff=0.0, width=50)
                    except Exception as exception:
                        print(exception)
    
                except Exception as exception:
                    print('\nRE-PLOT ERROR:', exception, '\n')
    
    # copy files
    if False:
        cwd = '/mnt/home/paterno1/Ben/otherCellTypes/'
        #cwd = os.path.dirname(__file__).replace('\\', '/') + '/'

        for workingDir in ['otherCellTypes_rec/', 'otherCellTypes_lig/']:
            for celltype in ['Endothelial', 'Macrophage', 'Pericyte', 'SMC', 'Fibroblast']:
                for species in ['Homo sapiens', 'Mus musculus']:
                    #for file in ['results All correlation euclidean ward.png',
                    #            'results All correlation euclidean ward.xlsx',
                    #            'all peaks Avg combo3avgs.png',
                    #            'All peaks Avg combo3avgs.xlsx',
                    #            'Avg combo3avgs_variant.xlsx']:    #'near frequency Avg combo3avgs.xlsx' dendrogram-heatmap-correlation-data.xlsx
                    for file in ['near frequency Avg combo3avgs.xlsx']:
                        source = cwd + workingDir + '%s/%s/' % (celltype, species)
                        destination = cwd + 'results 07 08 2021/cutoff_0.05/' + workingDir
                        #print(source)
                        #print(destination)
                        if not os.path.exists(destination):
                            os.makedirs(destination)

                        shutil.copyfile(source + file, destination + '%s %s. ' % (celltype, species) + file) #'bootstrap/All/'

                        if False:
                            if workingDir=='otherCellTypes_lig/' and file=='results All correlation euclidean ward.xlsx':
                            #if workingDir=='otherCellTypes_rec/' and file=='results All correlation euclidean ward.xlsx':
                                df = pd.read_excel(destination + '%s %s. ' % (celltype, species) + file, index_col=0, header=0)[['Markers']]
                                print('\n\n', workingDir, '\t', species, '\t', celltype, '\t', len(df))
                        
                                for gene in ['VEGFA']:
                                #for gene in ['KDR', 'FLT1']:
                                    print('\n%s:' % gene)
                        
                                    if gene in df.index.values:
                                        labels = df['Markers'].values
                                        values = np.abs(np.arange(len(df.index.values)) - np.where(df.index.values==gene)[0])
                        
                                        getNonParametricPValue(labels, values)
                                        getLogisticRegressionPValue(labels, values)
                                    else:
                                        print('  %s not in index' % gene)

