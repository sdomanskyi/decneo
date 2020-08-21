from common import *
from analysisPipeline import analyze, compareTwoCases

workingDir = 'otherCellTypes/'

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

    uCelltypes = ['Dendritic', 'Fibroblast', 'Erythroid', 'B', 'T', 'NK', 'Endothelial']
    dictCelltypes = {'Endothelial cells':                           'Endothelial', 
                'Fibroblasts':                                      'Fibroblast', 
                'Erythroid-like and erythroid precursor cells':     'Erythroid', 
                'T memory cells':                                   'T', 
                'B cells':                                          'B', 
                'T cells':                                          'T', 
                'NK cells':                                         'NK', 
                'Dendritic cells':                                  'Dendritic', 
                'Plasmacytoid dendritic cells':                     'Dendritic', 
                'Gamma delta T cells':                              'T', 
                'Endothelial cells (aorta)':                        'Endothelial', 
                'T regulatory cells':                               'T', 
                'B cells naive':                                    'B', 
                'T cells naive':                                    'T', 
                'T cytotoxic cells':                                'T'}

    # Load all selected cells IDs into one dataframe
    if False:
        from PanglaoDBannotation import getAnnotationsSummaryDf

        df_cell_type_annotations = getAnnotationsSummaryDf(MetadataDirName)[['Cell type annotation', 'Species']]
        df_cell_type_annotations = df_cell_type_annotations.loc[np.isin(df_cell_type_annotations['Cell type annotation'].values, list(dictCelltypes.keys()))]
        df_cell_type_annotations['Cell type annotation'] = df_cell_type_annotations['Cell type annotation'].replace(dictCelltypes)
        df_cell_type_annotations = df_cell_type_annotations.set_index(['Cell type annotation', 'Species'], append=True)
        print(df_cell_type_annotations)
    
        dfs = []
        for SRA, SRS, cluster, celltype, species in df_cell_type_annotations.index[:]:
            df_cell_to_cluster = pd.read_csv(os.path.join(MetadataDirName, 'PanglaoDB', 'data', 'sample_clusters', '%s%s.seurat_clusters.txt' % (SRA, '_' + SRS if SRS!='notused' else '')), delimiter=' ', index_col=0, header=None)[1]
            selectedCells = df_cell_to_cluster[df_cell_to_cluster == cluster].index.values
            dfs.append(pd.DataFrame(index=pd.MultiIndex.from_tuples([(species, '%s%s.sparse.RData.h5' % (SRA, '_' + SRS if SRS!='notused' else ''))], names=['species', 'file']), columns=np.unique(selectedCells), data=celltype).stack())

        dfs = pd.concat(dfs, axis=0, sort=False)
        dfs.to_hdf('data/selectedCellsInPanglaoDB.h5', key='df', mode='a', complevel=4, complib='zlib')
        print(dfs)

        if False:
            df = pd.read_hdf('data/selectedCellsInPanglaoDB.h5', key='df').to_frame()
            df.index.names = ['species', 'file', 'cell']
            df.columns = ['celltype']
            df = df.reset_index().set_index(['species', 'celltype'])
            df['file'] += '_' + df['cell']
            df = df['file']
            df = df.groupby(['species', 'celltype']).count().unstack('celltype')
            print(df)

    bigDataFile = processedDataDir + 'selectedCellsInPanglaoDB_processed.h5'

    # Prepare expression and DEG ranking.This step took 10 hrs
    if False:
        df_cells_id = pd.read_hdf('data/selectedCellsInPanglaoDB.h5', key='df')
        for species in ['Homo sapiens', 'Mus musculus']:
            df_cells_id_species = df_cells_id.xs(species, level='species')
            files = np.unique(df_cells_id_species.index.get_level_values('file'))
            print('Species:', species, '\tFiles:', len(files), '\tCells:', len(df_cells_id_species), flush=True)
            #np.savetxt('file_cell_id_%s.txt' % species, files, fmt='%s')

            for file in files[:]:
                try:
                    print('\nProcessing file:', file, flush=True)

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
                    df_expr = df_expr.astype(float)

                    # Scale and log-transform
                    df_expr /= df_expr.sum(axis=0) * 0.0001
                    df_expr = np.log2(df_expr.replace(0., np.min(df_expr.values[df_expr.values > 0.])))
                    df_expr -= np.min(df_expr.values)

                    # Remove nearly-constant and constant genes
                    df_expr = df_expr[np.std(df_expr, axis=1) / np.mean(np.std(df_expr.values)) > 0.01]

                    for celltype in np.unique(df_temp_cells.values):
                        celltypeCells = df_temp_cells[df_temp_cells == celltype].index
                        df_celltype = df_expr[celltypeCells]
                
                        # Rank and save differentially expressed genes
                        df_ttest = pd.DataFrame(index=df_celltype.index, columns=['statistic', 'pvalue'])
                        ttest = scipy.stats.ttest_ind(df_celltype.values, df_expr[df_expr.columns.difference(celltypeCells)].values, axis=1)
                        df_ttest['statistic'] = ttest[0]
                        df_ttest['pvalue'] = ttest[1]
                        df_ttest = df_ttest.sort_values('statistic', ascending=False).dropna()
                        df_ttest.to_hdf(saveFile, key='ttest/%s/%s/%s/' % (species, celltype, file), mode='a', complevel=4, complib='zlib')

                        # Save expression data
                        df_celltype = df_celltype[df_celltype.sum(axis=1) > 0]
                        df_celltype.to_hdf(bigDataFile, key='expr/%s/%s/%s/' % (species, celltype, file), mode='a', complevel=4, complib='zlib')

                        print('\t', celltype, df_ttest.shape, df_celltype.shape, flush=True)

                except Exception as exception:
                    print(exception, flush=True)

    # Aggregate files
    if False:
        from io import KeysOfStore

        keys = KeysOfStore(bigDataFile) # ['/ttest/Mus musculus/Endothelial/SRA203368_SRS866906.sparse.RData.h5', ...]
        skeys = np.array([np.array(key.strip('/').split('/')) for key in keys])

        if False:
            dfc = pd.DataFrame(skeys).set_index([0,1,2]).groupby([0,1,2]).count().unstack(2)
            dfc.to_excel('checkBigPerCelltype.xlsx')
            print(dfc)

            dfc = pd.DataFrame(skeys).set_index([0,1,3])
            dfc[:] = 1
            dfc.index.names = ['mode', 'species', 'file']
            dfc = dfc.loc[~dfc.index.duplicated(keep='first')]
            dfc = dfc.groupby(['mode', 'species']).count() #.unstack(3)
            dfc.to_excel('checkBigFiles.xlsx')
            print(dfc)

        dfc = pd.DataFrame(skeys).set_index([0,1,2])
        dfc.index.names = ['x', 'species', 'celltype']
        dfc.columns = ['file']

        for x in ['ttest', 'expr'][1:]:
            dfc_x = dfc.xs(x, level='x')
            for species in ['Homo sapiens', 'Mus musculus']:
                dfc_x_species = dfc_x.xs(species, level='species')
                for celltype in uCelltypes:
                    try:
                        tempFiles = dfc_x_species.xs(celltype).values.flatten().tolist()

                        if x == 'ttest':
                            genes = []
                            batches = []
                            for file in tempFiles:
                                df_temp = pd.read_hdf(bigDataFile, key='%s/%s/%s/%s/' % (x, species, celltype, file))
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
                            dfs = []
                            for file in tempFiles:
                                batch = file.split('.sparse.RData.h5')[0]
                                df_temp = pd.read_hdf(bigDataFile, key='%s/%s/%s/%s/' % (x, species, celltype, file))
                                df_temp = pd.concat([df_temp], keys=[batch], axis=1, sort=False)
                                df_temp.columns.names = ['batch', 'cell']
                                dfs.append(df_temp)

                            print('\nMerging expression data:', species, celltype, flush=True)
                            df_expr = pd.concat(dfs, axis=1, sort=False).fillna(0.)
                            
                            del dfs[:]

                            print('Saving expression data:', species, celltype, df_expr.shape, flush=True)
                            df_expr.to_hdf(processedDataDir + 'PanglaoDB_expresion_per_batch_%s_%s.h5' % (species, celltype), key='df', mode='a', complevel=4, complib='zlib')

                    except Exception as exception:
                        print('\nERROR:', exception, '\n')

    # Analyze per celltype
    if True:
        try:
            batchID = eval(sys.argv[1]) # 0, 1, 2
            
            if batchID == 0:
                uCelltypes = ['Fibroblast', 'Erythroid']
            elif batchID == 1:
                uCelltypes = ['B', 'T']
            elif batchID == 2:
                uCelltypes = ['NK', 'Endothelial']

            print('batchID:', batchID, bootstrapExperiments, flush=True)
        except:
            pass

        for celltype in uCelltypes:
            try:
                comparisonName = workingDir + '%s/%s/comparison' % (celltype, 'Mus musculus')

                if True:
                    for species in ['Homo sapiens', 'Mus musculus']:
                        print('Analyzing all cells of %s, %s:' % (species, celltype))

                        df_expr = pd.read_hdf(processedDataDir + 'PanglaoDB_expresion_per_batch_%s_%s.h5' % (species, celltype), key='df')

                        analyze(df_expr, receptorsListHugo_2555, gECs, gECi, 'correlation',
                                suffix='%s_%s' % (species, celltype), saveDir=workingDir + '%s/%s/' % (celltype, species),
                                toggleCalculateMajorMetric=True, toggleExportFigureData=True, toggleCalculateMeasures=True,
                                toggleAdjustText=False, noPlot=True, panels=[])

                    compareTwoCases(workingDir + '%s/%s/' % (celltype, 'Homo sapiens'), 
                                    workingDir + '%s/%s/' % (celltype, 'Mus musculus'), 
                                    saveName=comparisonName)
    
                for species in ['Homo sapiens', 'Mus musculus']:
                    print('Re-analyzing for of %s, %s:' % (species, celltype))

                    additionalData = externalPanelsData.copy()
                    additionalData.update({
                        'diffExpressedGenes': pd.read_hdf(processedDataDir + 'PanglaoDB_ttest_ranks_per_batch_%s_%s.h5' % (species, celltype), key='df'),
                        'conservedGenes': pd.read_excel(comparisonName + '.xlsx', index_col=1, header=0)['Inter-measures.T50_common_count']})

                    analyze(None, receptorsListHugo_2555, gECs, gECi, 'correlation',
                            suffix='%s_%s' % (species, celltype), saveDir=workingDir + '%s/%s/' % (celltype, species), toggleIncludeHeatmap=True,
                            toggleCalculateMajorMetric=False, toggleExportFigureData=True, toggleCalculateMeasures=False,
                            externalPanelsData=additionalData) 
            
            except Exception as exception:
                print('\nERROR:', exception, '\n')

