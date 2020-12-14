from decneo.commonFunctions import *
from decneo.analysisPipeline import Analysis, process

def preparePerSampleDEgenes(wdir):

    allCellsIndexAnnotation = wdir + 'AlonaCellsInPanglaoDB.h5'
    ECIndexAnnotation = wdir + 'AlonaECinPanglaoDB.h5'
    nonECIndexAnnotation = wdir + 'AlonaNonECinPanglaoDB.h5'
    fileNameEC = wdir + 'PanglaoDB_Alona_EC.h5'
    fileNameOther = wdir + 'PanglaoDB_Alona_nonEC.h5'
    fileNameTtest = wdir + 'ttest_Alona_PanglaoEC.h5'
    fileNameRanks = wdir + 'ranks_Alona_PanglaoEC.h5'
    fileNameDfRanks = wdir + 'PanglaoDB_Alona_ttest_ranks_per_batch.h5'

    # Prepare index of all cells
    if False:
        df_cell_type_annotations = getPanglaoDBAnnotationsSummaryDf(MetadataDirName)[['Cell type annotation', 'Species']]
        df_cell_type_annotations = df_cell_type_annotations.set_index(['Cell type annotation', 'Species'], append=True)
        df_cell_type_annotations.index.names = ['SRA', 'SRS', 'cluster', 'celltype', 'species']
        print(df_cell_type_annotations.index.names)
        print(df_cell_type_annotations.index[0])
    
        dfs = []
        for i, (SRA, SRS, cluster, celltype, species) in enumerate(df_cell_type_annotations.index):
            print(i, 'of', len(df_cell_type_annotations), SRA, SRS, cluster, celltype, species)
            se_cells = pd.read_csv(os.path.join(MetadataDirName, 'PanglaoDB', 'data', 'sample_clusters', '%s%s.seurat_clusters.txt' % (SRA, '_' + SRS if SRS!='notused' else '')), delimiter=' ', index_col=0, header=None)[1]
            se_cells = se_cells.to_frame().reset_index().set_index(1)
            se_cells.index.names = ['cluster']
            se_cells = se_cells.xs(cluster, drop_level=False)

            se_cells = pd.concat([se_cells], keys=[species], names=['species'], sort=False)
            se_cells = pd.concat([se_cells], keys=[celltype], names=['celltype'], sort=False)
            se_cells = pd.concat([se_cells], keys=[SRS], names=['SRS'], sort=False)
            se_cells = pd.concat([se_cells], keys=[SRA], names=['SRA'], sort=False)
            se_cells = se_cells.reorder_levels(['SRA', 'SRS', 'cluster', 'celltype', 'species'])
            se_cells.columns = ['cell']

            se_cells['file'] = '%s%s.sparse.RData.h5' % (SRA, '_' + SRS if SRS!='notused' else '')
            
            dfs.append(se_cells)
            
        dfs = pd.concat(dfs, axis=0, sort=False)
        dfs.to_hdf(allCellsIndexAnnotation, key='df', mode='a', complevel=4, complib='zlib')
        print(dfs)
        
    # Prepare index of EC and non-EC   
    if False:
        df_all_cells = pd.read_hdf(allCellsIndexAnnotation, key='df')

        df_ECindex = df_all_cells.loc[np.isin(df_all_cells.index.get_level_values('celltype'), ['Endothelial cells', 'Endothelial cells (aorta)'])]
        df_ECindex.to_hdf(ECIndexAnnotation, key='df', mode='a', complevel=4, complib='zlib')
        print(df_ECindex)

        df_nonECindex = df_all_cells.loc[~np.isin(df_all_cells.index.get_level_values('celltype'), ['Endothelial cells', 'Endothelial cells (aorta)'])]
        print(df_nonECindex)

        dfs = []
        for i, (SRA, SRS) in enumerate(pd.MultiIndex.from_arrays([df_ECindex.index.get_level_values('SRA'), df_ECindex.index.get_level_values('SRS')]).to_frame().drop_duplicates().index.values):
            print(i, SRA, SRS)
            df_temp = df_nonECindex.xs(key=SRA, level='SRA', drop_level=False).xs(key=SRS, level='SRS', drop_level=False)
            df_temp = df_temp.sample(min(1000, df_temp.shape[0]), replace=False)
            dfs.append(df_temp)

        dfs = pd.concat(dfs, axis=0, sort=False)
        dfs.to_hdf(nonECIndexAnnotation, key='df', mode='a', complevel=4, complib='zlib')
        print(dfs)
     
    # Collect expression of EC and non-EC and determine ranks of DE genes
    if True:
        df_ECindex = pd.read_hdf(ECIndexAnnotation, key='df')
        df_nonECindex = pd.read_hdf(nonECIndexAnnotation, key='df')

        print(df_ECindex)
        print(df_nonECindex)

        files = np.unique(df_ECindex['file'])

        for i, file in enumerate(files):
            try:
                batch = file.split('.')[0]

                if KeyInStore(batch, fileNameRanks):
                    continue

                temp = batch.split('_')
                if len(temp)==2:
                    SRA, SRS = temp
                else:
                    SRA, SRS = temp[0], 'notused'

                print('Processing:', batch, '\t', i, 'of', len(files), flush=True)
                df_expr = pd.read_hdf(RDataDirName + file, key='df')
                print(df_expr.shape)

                cellsEC = df_ECindex.xs(key=SRA, level='SRA', drop_level=False).xs(key=SRS, level='SRS', drop_level=False)['cell'].values

                try:
                    cellsNonEC = df_nonECindex.xs(key=SRA, level='SRA', drop_level=False).xs(key=SRS, level='SRS', drop_level=False)['cell'].values
                except:
                    cellsNonEC = []
                
                if len(cellsEC) >= 10:
                    df_expr = pd.concat([df_expr[cellsEC], df_expr[cellsNonEC]], axis=1, sort=False).fillna(0.)
                    df_expr.index = pd.Series(df_expr.index.values).replace(Mouse_to_Human_HUGO_conversion).values
                    df_expr = df_expr.loc[~df_expr.index.duplicated(keep='first')]
                    df_expr = df_expr.T.loc[~df_expr.T.index.duplicated(keep='first')].T
                    df_expr = df_expr.astype(float).fillna(0.)
                    df_expr = df_expr[df_expr.columns[df_expr.sum(axis=0) > 0.]]
                    df_expr /= df_expr.sum(axis=0) * 0.0001
                    df_expr = np.log2(df_expr.replace(0., np.min(df_expr.values[df_expr.values > 0.])))
                    df_expr -= np.min(df_expr.values)

                    df_EC, df_other = df_expr[cellsEC], df_expr[cellsNonEC]

                    df_EC = pd.concat([df_EC], keys=[batch], names=['batch'], axis=1, sort=True)
                    df_EC.columns.names = ['batch', 'cell']

                    df_other = pd.concat([df_other], keys=[batch], names=['batch'], axis=1, sort=True)
                    df_other.columns.names = ['batch', 'cell']

                    df_EC.to_hdf(fileNameEC, key=batch, mode='a', complevel=4, complib='zlib')
                    df_other.to_hdf(fileNameOther, key=batch, mode='a', complevel=4, complib='zlib')
                
                    df_ttest = pd.DataFrame(index=df_EC.index, columns=['statistic', 'pvalue'])
                    ttest = scipy.stats.ttest_ind(df_EC.values, df_other.values, axis=1)
                    df_ttest['statistic'] = ttest[0]
                    df_ttest['pvalue'] = ttest[1]
                    df_ttest = df_ttest.sort_values('statistic', ascending=False).dropna()
                    df_ttest.to_hdf(fileNameTtest, key=batch, mode='a', complevel=4, complib='zlib')

                    se_ranks = pd.Series(data=df_ttest.loc[df_ttest['pvalue'] <= 10**-3]['statistic'].index.values)
                    se_ranks.to_hdf(fileNameRanks, key=batch, mode='a', complevel=4, complib='zlib')

                    print(df_EC.shape[0], df_EC.shape[1], df_other.shape[1], df_ttest.shape, se_ranks.shape, flush=True)

            except Exception as exception:
                print('Error while processing file:\n', exception)

    # Prepare df_ranks
    if True:
        allBatches = pd.read_hdf(ECIndexAnnotation, key='df')['file'].droplevel(['cluster', 'celltype', 'SRA', 'SRS']).str.split('.', expand=True)[0].drop_duplicates()
        
        for species in ['Mus musculus', 'Homo sapiens']:
            thisSpeciesBatches = allBatches.xs('Homo sapiens').values if species == 'Homo sapiens' else allBatches.xs('Mus musculus').values

            df_ranks = []
            for batch in np.intersect1d(np.array([key[1:] for key in KeysOfStore(fileNameRanks)]), thisSpeciesBatches):
                se_ranks = pd.read_hdf(fileNameRanks, key=batch)
                df_ranks.append(pd.Series(index=se_ranks.values, data=se_ranks.index, name=batch))

            df_ranks = pd.concat(df_ranks, axis=1, sort=False).add(1.).fillna(-1).astype(int)
            df_ranks.to_hdf(fileNameDfRanks, key=species, mode='a', complevel=4, complib='zlib')
            print(df_ranks)

    return

def prepareInput(wdir, fileName, species):

    df_ranks = pd.read_hdf(wdir + 'PanglaoDB_Alona_ttest_ranks_per_batch.h5', key=species)
    df_ranks.to_hdf(wdir + fileName, key='df_ranks', mode='a', complevel=4, complib='zlib')
    print(df_ranks)

    df = pd.concat([pd.read_hdf(wdir + 'PanglaoDB_Alona_EC.h5', key=batch) for batch in df_ranks.columns], axis=1, sort=False).fillna(0.)
    df.to_hdf(wdir + fileName, key='df', mode='a', complevel=4, complib='zlib')
    print(df)

    return

if __name__ == '__main__':

    if platform.system()=="Windows":
        wdir = 'd:/Projects/A_Endothelial/VS/Endothelial/results/PanglaoDB_byAlona/' 
    else:
        wdir = '/mnt/research/piermarolab/Sergii/PanglaoDB_byAlona/'
    
    if False:
        preparePerSampleDEgenes()
        prepareInput(wdir, 'PanglaoDB_byDCS_human/data.h5', 'Homo sapiens')
        prepareInput(wdir, 'PanglaoDB_byDCS_mouse/data.h5', 'Mus musculus')

    #anHuman, anMouse = process(*(None, None), *(None, None),
    #                            wdir + 'PanglaoDB_byDCS_human/', wdir + 'PanglaoDB_byDCS_mouse', 
    #                            nCPUs=4 if platform.system()=="Windows" else 20, parallelBootstrap=False,
    #                            PCNpath=os.path.join(os.path.dirname(__file__), 'data'), exprCutoff1=0.05, 
    #                            genesOfInterest=receptorsListHugo_2555, knownRegulators=gEC22, perEachOtherCase=True)
    #
    #anHuman.reanalyzeMain(togglePublicationFigure=True, toggleIncludeHeatmap=False, markersLabelsRepelForce=1.25)
    #anMouse.reanalyzeMain(togglePublicationFigure=True, toggleIncludeHeatmap=False, markersLabelsRepelForce=1.5)
    
    if True:
        kwargs = dict(genesOfInterest=receptorsListHugo_2555, knownRegulators=gEC22, 
                      panels=[
                          'combo3avgs', 
                          'combo4avgs', 
                          'fraction', 
                          'binomial', 
                          'top50',
                          'markers', 
                          ], perEachOtherCase=True, PCNpath=os.path.join(os.path.dirname(__file__), 'data'))

        anHuman = Analysis(workingDir=wdir + 'PanglaoDB_byDCS_human/', otherCaseDir=wdir + 'PanglaoDB_byDCS_mouse/', **kwargs)
        anHuman.reanalyzeMain(togglePublicationFigure=True, toggleIncludeHeatmap=False, markersLabelsRepelForce=1.25, includeClusterNumber=False)

        anMouse = Analysis(workingDir=wdir + 'PanglaoDB_byDCS_mouse/', otherCaseDir=wdir + 'PanglaoDB_byDCS_human/', **kwargs)
        anMouse.reanalyzeMain(togglePublicationFigure=True, toggleIncludeHeatmap=False, markersLabelsRepelForce=1.5, includeClusterNumber=False)