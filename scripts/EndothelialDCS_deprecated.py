from scRegulation.commonFunctions import *

receptors = np.loadtxt('geneLists/receptors_Human_HUGO.txt', dtype=str)
ligands = np.loadtxt('geneLists/ligands_Human_HUGO.txt', dtype=str)
Mouse_to_Human_HUGO_conversion = read('geneLists/Mouse_to_Human_HUGO.gz', jsonFormat=True)
tip = pd.Series(np.loadtxt('geneLists/tip.txt', dtype=str)).replace(Mouse_to_Human_HUGO_conversion).values.tolist()
stalk = pd.Series(np.loadtxt('geneLists/stalk.txt', dtype=str)).replace(Mouse_to_Human_HUGO_conversion).values.tolist()
df_TF = pd.read_excel('geneLists/TranscriptionFactors_DatabaseExtract_v_1.01.xlsx', header=0)[['HGNC symbol', 'Is TF?']]
listTF = df_TF['HGNC symbol'][df_TF['Is TF?'] == 'Yes'].values.tolist()
listUnlikelyTF = df_TF['HGNC symbol'][df_TF['Is TF?'] == 'No'].values.tolist()

if False:
    HGNC_receptors = np.loadtxt('HGNC_receptors.txt', dtype=str)
    HGNC_receptor_ligands = pd.read_excel('HGNC receptor ligands SD.xlsx', index_col=0, header=0).index.drop_duplicates().values.tolist()
    recAlmen = pd.read_excel('Almen 12915_2009_258_MOESM1_ESM copy.xls', sheet_name='Sheet3', header=0)['Rec Almen'].dropna()
    recAlmen = recAlmen.str.split(';', expand=True).stack().droplevel(1).drop_duplicates().values.tolist()
    recRamilow = pd.read_excel('Almen 12915_2009_258_MOESM1_ESM copy.xls', 
                            sheet_name='Sheet3', header=0)['Rec Ramilow'].dropna().drop_duplicates().values.tolist()

    recList = np.unique(receptors.tolist() + HGNC_receptors.tolist() + recAlmen + recRamilow).tolist() # 2634
    DCS = DigitalCellSorter.DigitalCellSorter()
    recListHugo = DCS.gnc.Convert(recList, sourceType='alias', targetType='hugo', returnUnknownString=False)
    pd.Series(np.unique(recListHugo)).to_excel('receptorsListHugo.xlsx')
    exit()

if False:
    allHGNCgenes = pd.read_csv('all HGNC genes with groups.txt', index_col=None, header=0, delimiter='\t')[['Approved symbol', 'Status']].set_index('Approved symbol')['Status']
    allHGNCgenes = allHGNCgenes[allHGNCgenes == 'Approved'].index.values.tolist()

    allLists = ['receptors', 'ligands', 'HGNC_receptors', 'HGNC_receptor_ligands', 'recAlmen', 'recRamilow', 'tip', 'stalk', 'listTF', 'listUnlikelyTF', 'allHGNCgenes']
    df = pd.DataFrame(index=allLists, columns=allLists)
    for geneList in allLists:
        geneList = set(np.unique(eval(geneList)))
        print(len(geneList))

    for col in df.columns:
        for row in df.index:
            set1 = set(np.unique(eval(row)))
            set2 = set(np.unique(eval(col)))
            df.loc[row, col] = len(set1.intersection(set2))

    df.to_excel('geneLists comparison SD.xlsx')
    print(df)


if __name__ == '__main__':

    np.random.seed(0)

    # Combine samples to one DataFrame
    if False:
        dataFile = 'DCS output/PanglaoDBendothelialAllv0.h5'

        keys = KeysOfStore(dataFile)
        print('Reading %s keys of %s' % (len(keys), dataFile), flush=True)

        dfs = []
        counts = 0
        for ikey, key in enumerate(keys):
            key = key[1:]

            print('\tReading:', ikey, 'of', len(keys), key, 'Count:', counts, flush=True)
            df_temp = pd.read_hdf(dataFile, key=key).astype('float32')

            if df_temp.shape[1] >= 20:
                dfs.append(df_temp)

                print(df_temp.shape, df_temp.dtypes[0], flush=True)
                counts += df_temp.shape[1]

        print('Total count:', counts, flush=True)

        print('Combining samples', flush=True)
        dfs = pd.concat(dfs, axis=1, sort=False).fillna(0.)
        print(dfs, dfs.dtypes[0], flush=True)

        dfs = dfs.loc[(dfs > 0).sum(axis=1) >= 100]
        print(dfs, dfs.dtypes[0], flush=True)

        dfs.to_hdf(dataFile[:-3] + '_combined.h5', key='df', mode='a', complevel=4, complib='zlib')

    # Analyze Human and Mouse samples (separately)
    if True:
        dataFile = 'DCS output/PanglaoDBendothelialAllv0_combined.h5'

        for species in ['Homo sapiens', 'Mus musculus']:
            DCS = DigitalCellSorter.DigitalCellSorter(dataName='PanglaoDB_EC_DCS', saveDir='DCS output/PanglaoDB_EC_DCS/%s/'%species, verbose=2)

            # Get PnaglaoDB metadata (df_anno)
            if False:
                fields = ['Tissue origin of the sample', 'Species', 'Number of raw cells', 'Fraction of cells passed QC', 'Tissue origin of the sample', 'Sequencing instrument', 'scRNA-seq protocol', 'Median number of expressed genes per cell', 'Is the sample from primary adult tissue?', 'Is the sample from a tumor? (1 true otherwise false)', 'Is the sample from a cell line?']
                df_anno = getAnnotationsSummaryDf(MetadataDirName).droplevel('Cluster index', axis=0)[fields]
                df_anno.columns = ['Tissue', 'Species', 'Sequenced cells', 'QC-passed fraction', 'Tissue detailed', 'Sequencing instrument', 'scRNA-seq protocol', 'medianOfGenesPerCell', 'isPrimaryAdult', 'isTumor', 'isCellLine']
                df_anno = df_anno.loc[~df_anno.index.duplicated(keep='first')]
                df_anno.index = ['%s%s' % (item[0], '_' + item[1] if item[1]!='notused' else '') for item in df_anno.index]

                tissueDict = np.loadtxt('Tissue_dict_PanglaoDB_bySD.txt', delimiter='\t', dtype=str)
                tissueDict = dict(zip(tissueDict.T[1], tissueDict.T[0]))
                df_anno['Tissue'] = df_anno['Tissue'].replace(tissueDict)
                print(df_anno)

            # Step 1. Un-cut EC clustering and analysis
            if False:
                # Project, cluster, plot with DCS yet "un-cut" Endothelial cells
                if False:
                    df = pd.read_hdf(dataFile, key='df')
                    df = df[df.columns[np.isin(df.columns.get_level_values('batch').values, df_anno[df_anno['Species'] == species].index)]]
                    df = df.T[df.sum(axis=0) > 0].T
                    df = df[df.sum(axis=1) > 0]
                    print(df)

                    DCS.prepare(df.astype(float))

                    del df

                    if False:
                        DCS.project()

                    if False:
                        DCS.nClusters = 15
                        DCS.cluster()
                        DCS.makeProjectionPlotByClusters(suffix='by Clusters')

                    print('Making gene plot(s)')
                    df_projection = pd.read_hdf(DCS.fileHDFpath, key='df_projection')

                    valuesEndo = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Endothelial').reindex(df_projection.columns)
                    valuesEpi = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Epithelial').reindex(df_projection.columns)
                    valuesFibro = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Fibroblasts').reindex(df_projection.columns)
                    vEndo = ((valuesEndo.fillna(0.).values >= 0.2) & (valuesEpi.fillna(0.).values == 0.) & (valuesFibro.fillna(0.).values == 0.)).astype(bool)
                    print(len(vEndo), vEndo.sum(), vEndo.sum()/len(vEndo))
                    df_projection = df_projection[df_projection.columns[vEndo]]

                    DCS.geneListFileName = 'geneLists/fibro_endo_epi_v2_Human.xlsx'
                    markers = DCS.readMarkerFile().index.values.tolist()

                    print('Markers:', len(markers))

                    for gene in markers + ['FLT4', 'VEGFA', 'VEGFB', 'VEGFC', 'VEGFD']:
                        try:
                            expr = DCS.df_expr.xs(key=DCS.getHugoName(gene)).loc[df_projection.columns].replace(0., np.nan).values
                            DCS.makeProjectionPlot(df_projection.values, reduce(expr), legend=False, labels=False, colorbar=True, suffix='by %s' % gene)
                        except:
                            pass

                    df_projection = pd.read_hdf(DCS.fileHDFpath, key='df_projection')

                    valuesEndo = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Endothelial').reindex(df_projection.columns)
                    valuesEpi = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Epithelial').reindex(df_projection.columns)
                    valuesFibro = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Fibroblasts').reindex(df_projection.columns)
                    vEndo = ((valuesEndo.fillna(0.).values >= 0.2) & (valuesEpi.fillna(0.).values == 0.) & (valuesFibro.fillna(0.).values == 0.)).astype(bool)
                    print(len(vEndo), vEndo.sum(), vEndo.sum()/len(vEndo))
                    df_projection = df_projection[df_projection.columns[vEndo]]

                    df_anno = df_anno.reindex(df_projection.columns.get_level_values('batch'), axis=0)

                    for col in ['Tissue', 'Species', 'Tissue detailed', 'Sequencing instrument', 'scRNA-seq protocol']:
                        DCS.makeProjectionPlot(df_projection.values, df_anno[col].values, legend=True, labels=True, suffix='by %s' % col)

                    for col in ['Sequenced cells', 'QC-passed fraction']:
                        DCS.makeProjectionPlot(df_projection.values, reduce(df_anno[col].values), legend=False, labels=False, colorbar=True, suffix='by %s' % col)

                # Projection plots and Sankey
                if False:
                    # Prepare PanglaoDB celltype annotations for comparison with DCS annotations
                    if False:
                        annotationPath = os.path.join(MetadataDirName, 'PanglaoDB', 'data', 'cell_type_annotations.txt')
                        df_cell_type_annotations = pd.read_csv(annotationPath, index_col=[0, 1, 2], header=None)[3]
                        df_cell_type_annotations.index.names = ['SRA', 'SRS', 'cluster']
                        df_cell_type_annotations.name = 'celltype'
                        df_cell_type_annotations = df_cell_type_annotations.reset_index()
                        df_cell_type_annotations.index = [i + ('_' + j if j != 'notused' else '') for i,j in zip(df_cell_type_annotations['SRA'].values, df_cell_type_annotations['SRS'].values)]
                        df_cell_type_annotations = df_cell_type_annotations.drop(['SRA', 'SRS'], axis=1)
                        df_cell_type_annotations.index.name = 'batch'
                        df_cell_type_annotations = df_cell_type_annotations.reset_index().set_index(['batch', 'cluster'])
                        print(df_cell_type_annotations)

                        cells = pd.read_hdf(DCS.fileHDFpath, key='df_projection').columns
            
                        seAnno = []
                        for ibatch, batch in enumerate(cells.levels[0].values[:]):
                            print(batch)
                            p = batch.split('_')
                            SRA, SRS = p[0], 'notused' if len(p)==1 else p[1]

                            batchCells = pd.read_csv(os.path.join(MetadataDirName, 'PanglaoDB', 'data', 'sample_clusters', '%s%s.seurat_clusters.txt' % (SRA, '_' + SRS if SRS!='notused' else '')), delimiter=' ', index_col=0, header=None)[1]
                            batchCells = batchCells.loc[~batchCells.index.duplicated(keep='first')]
                            batchCells = pd.concat([batchCells], keys=[batch])
                            batchCells = batchCells.reset_index()
                            batchCells.columns = ['batch', 'cell', 'cluster']
                            batchCells = batchCells.set_index(['batch', 'cluster'])
                            batchCells['celltype'] = df_cell_type_annotations.reindex(batchCells.index).fillna('Missing').values
                            batchCells = batchCells.reset_index().set_index(['batch', 'cell'])['celltype']

                            seAnno.append(batchCells)

                        seAnno = pd.concat(seAnno, axis=0, sort=False)
                        seAnno.to_hdf('data/seAnnotaionPanglaoDBsubsetEC.h5', key='seAnno/%s' % species, mode='a')
                        print(seAnno)
                    else:
                        seAnno = pd.read_hdf('data/seAnnotaionPanglaoDBsubsetEC.h5', key='seAnno/%s' % species)
                        print(seAnno)

                    df_projection = pd.read_hdf(DCS.fileHDFpath, key='df_projection')

                    valuesEndo = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Endothelial').reindex(df_projection.columns)
                    valuesEpi = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Epithelial').reindex(df_projection.columns)
                    valuesFibro = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Fibroblasts').reindex(df_projection.columns)
                    vEndo = ((valuesEndo.fillna(0.).values >= 0.2) & (valuesEpi.fillna(0.).values == 0.) & (valuesFibro.fillna(0.).values == 0.)).astype(bool)
                    print(len(vEndo), vEndo.sum(), vEndo.sum()/len(vEndo))
                    df_projection = df_projection[df_projection.columns[vEndo]]

                    seTissue = pd.Series(index=df_projection.columns, data=df_anno['Tissue'].reindex(df_projection.columns.get_level_values('batch')).fillna(0.).values)
                    seTissueDetailed = pd.Series(index=df_projection.columns, data=df_anno['Tissue detailed'].reindex(df_projection.columns.get_level_values('batch')).fillna(0.).values)
                    seBatch = pd.Series(index=df_projection.columns, data=df_projection.columns.get_level_values('batch').values)
                    sePanglaoCelltype = pd.Series(index=df_projection.columns, data=seAnno.reindex(df_projection.columns).fillna('Missing').values)
                    dfCounts = DCS.getCountsDataframe(seBatch, seTissue)
                    DCS.makeSankeyDiagram(dfCounts, attemptSavingHTML=True)

                    DCS.makeProjectionPlot(df_projection.values, vAnno, legend=True, labels=True, suffix='by is PnaglaoDB_Anno')

                # Projection by DCS celltype assignment confidence on "un-cut" version
                if False:
                    df_projection = pd.read_hdf(DCS.fileHDFpath, key='df_projection')

                    valuesEndo = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Endothelial').reindex(df_projection.columns)
                    valuesEpi = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Epithelial').reindex(df_projection.columns)
                    valuesFibro = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Fibroblasts').reindex(df_projection.columns)

                    vEndo = ((valuesEndo.fillna(0.).values >= 0.2) & (valuesEpi.fillna(0.).values == 0.) & (valuesFibro.fillna(0.).values == 0.)).astype(bool)
                    print(len(vEndo), vEndo.sum(), vEndo.sum()/len(vEndo))

                    DCS.makeProjectionPlot(df_projection.values, vEndo, legend=True, labels=True, suffix='by is Endothelial')
                    DCS.makeProjectionPlot(df_projection.values, valuesEndo.values, legend=False, labels=False, colorbar=True, suffix='by Endothelial confidence')
                    DCS.makeProjectionPlot(df_projection.values, valuesEpi.values, legend=False, labels=False, colorbar=True, suffix='by Epithelial confidence')
                    DCS.makeProjectionPlot(df_projection.values, valuesFibro.values, legend=False, labels=False, colorbar=True, suffix='by Fibroblasts confidence')
        
            # Step 2. Cut EC clustering and analysis
            if False:
                # Project and cluster with DCS (making "cut" filtered dataset of Endothelial cells)
                if False:
                    df_projection = pd.read_hdf(DCS.fileHDFpath, key='df_projection')
                    valuesEndo = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Endothelial').reindex(df_projection.columns)
                    valuesEpi = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Epithelial').reindex(df_projection.columns)
                    valuesFibro = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Fibroblasts').reindex(df_projection.columns)
                    cutoff = 0.75
                    vEndo = ((valuesEndo.fillna(0.).values >= cutoff) & (valuesEpi.fillna(0.).values == 0.) & (valuesFibro.fillna(0.).values == 0.)).astype(bool)
                    print(len(vEndo), vEndo.sum(), vEndo.sum()/len(vEndo))
                    df_projection = df_projection[df_projection.columns[vEndo]]
                    print(df_projection.shape)

                    DCS.saveDir = 'DCS output/PanglaoDB_EC_DCS/%s/cut/'%species

                    df = pd.read_hdf(dataFile, key='df')[df_projection.columns]
                    df = df.T[df.sum(axis=0) > 0].T
                    df = df[df.sum(axis=1) > 0]
                    print(df)

                    DCS.prepare(df.astype(float))
                    del df

                    if True:
                        DCS.project()

                    if True:
                        DCS.saveDir = 'DCS output/PanglaoDB_EC_DCS/%s/cut/'%species
                        DCS.nClusters = 15
                        DCS.cluster()
                        DCS.makeProjectionPlotByClusters(suffix='by Clusters')

                # Sankey and projection plots of "cut" version
                if False:
                    df_projection = pd.read_hdf(DCS.fileHDFpath, key='df_projection')

                    valuesEndo = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Endothelial').reindex(df_projection.columns)
                    valuesEpi = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Epithelial').reindex(df_projection.columns)
                    valuesFibro = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Fibroblasts').reindex(df_projection.columns)
                    cutoff = 0.75
                    vEndo = ((valuesEndo.fillna(0.).values >= cutoff) & (valuesEpi.fillna(0.).values == 0.) & (valuesFibro.fillna(0.).values == 0.)).astype(bool)
                    print(len(vEndo), vEndo.sum(), vEndo.sum()/len(vEndo))
                    df_projection = df_projection[df_projection.columns[vEndo]]

                    print(df_projection.shape)

                    ###
                    df_anno = df_anno.reindex(df_projection.columns.get_level_values('batch'), axis=0)

                    for col in ['Tissue', 'Species', 'Tissue detailed', 'Sequencing instrument', 'scRNA-seq protocol']:
                        DCS.makeProjectionPlot(df_projection.values, df_anno[col].values, legend=True, labels=True, suffix='by %s cut' % col)

                    for col in ['Sequenced cells', 'QC-passed fraction']:
                        DCS.makeProjectionPlot(df_projection.values, reduce(df_anno[col].values), legend=False, labels=False, colorbar=True, suffix='by %s cut' % col)

                    ###
                    seAnno = pd.read_hdf('data/seAnnotaionPanglaoDBsubsetEC.h5', key='seAnno/%s' % species)
                    seTissue = pd.Series(index=df_projection.columns, data=df_anno['Tissue'].reindex(df_projection.columns.get_level_values('batch')).fillna(0.).values)
                    seTissueDetailed = pd.Series(index=df_projection.columns, data=df_anno['Tissue detailed'].reindex(df_projection.columns.get_level_values('batch')).fillna(0.).values)
                    seBatch = pd.Series(index=df_projection.columns, data=df_projection.columns.get_level_values('batch').values)
                    sePanglaoCelltype = pd.Series(index=df_projection.columns, data=seAnno.reindex(df_projection.columns).fillna('Missing').values)
                    DCS.makeSankeyDiagram(DCS.getCountsDataframe(seBatch, seTissue), attemptSavingHTML=True, nameAppend=' Sankey batch-tissue')
                    DCS.makeSankeyDiagram(DCS.getCountsDataframe(seTissueDetailed, seTissue), attemptSavingHTML=True, nameAppend=' Sankey tissueDet-tissue')
                    DCS.makeSankeyDiagram(DCS.getCountsDataframe(seTissue, sePanglaoCelltype), attemptSavingHTML=True, nameAppend=' Sankey tissue-celltype')
                    DCS.makeSankeyDiagram(DCS.getCountsDataframe(seTissueDetailed, sePanglaoCelltype), attemptSavingHTML=True, nameAppend=' Sankey tissueDet-celltype')

                    ###
                    print('Plotting')
                    DCS.makeProjectionPlot(df_projection.values, reduce(valuesEpi.reindex(df_projection.columns)), legend=False, labels=False, colorbar=True, suffix='by is PnaglaoDB_confEpi_cut%s'%cutoff)
                    DCS.makeProjectionPlot(df_projection.values, reduce(valuesFibro.reindex(df_projection.columns)), legend=False, labels=False, colorbar=True, suffix='by is PnaglaoDB_confFibro_cut%s'%cutoff)
                    DCS.makeProjectionPlot(df_projection.values, reduce(valuesEndo.reindex(df_projection.columns)), legend=False, labels=False, colorbar=True, suffix='by is PnaglaoDB_confEndo_cut%s'%cutoff)
  
            # Step 3. Post-clustering plots and analysis of "cut" version
            if True:
                DCS = DigitalCellSorter.DigitalCellSorter(dataName='PanglaoDB_EC_DCS', saveDir='DCS output/PanglaoDB_EC_DCS/%s/cut/'%species, verbose=2)

                # Export for Michelle
                if False:
                    def getDf(species):

                        DCS = DigitalCellSorter.DigitalCellSorter(dataName='PanglaoDB_EC_DCS', saveDir='DCS output/PanglaoDB_EC_DCS/%s/cut/'%species, verbose=2)

                        df_clusters = pd.read_hdf(DCS.fileHDFpath, key='df_clusters')
                        df_clusters.index = df_clusters.index.get_level_values(0).values + df_clusters.index.get_level_values(1).values
                        df_clusters.index.name = None

                        pref = 'mouse' if species=='Mus musculus' else 'human'
                        df_clusters['cluster'] = pref + '-' + df_clusters['cluster'].str.replace('.0.','.').str.strip('.')
                        df_clusters = df_clusters.sort_values(by='cluster')
                        print(df_clusters)

                        return df_clusters

                    df_clusters = getDf('Homo sapiens').append(getDf('Mus musculus'))
                    print(df_clusters)
                    df_clusters.to_excel('Endothelial_DCS_clusters_mouse_and_human.xlsx')

                df_projection = pd.read_hdf(DCS.fileHDFpath, key='df_projection')

                #pd.read_hdf(DCS.fileHDFpath, key='df_clusters').to_hdf('DCS_EC_cut_df_clusters.h5', key=species)

                # Plots
                if False:
                    ### Gene heatmap plot
                    if False:
                        DCS.geneListFileName = 'geneLists/fibro_endo_epi_v2_Human.xlsx'
                        markers = DCS.readMarkerFile().index.values.tolist()
                        print('Markers:', len(markers))
                        selGenesHUGO = [DCS.getHugoName(gene) for gene in markers + ['FLT4', 'VEGFA', 'VEGFB', 'VEGFC', 'VEGFD'] + tip + stalk]

                        seCluster = pd.read_hdf(DCS.fileHDFpath, key='df_clusters').reindex(df_projection.columns)['cluster']

                        print('Reading data')
                        if False:
                            df = pd.read_hdf(dataFile, key='df')[df_projection.columns]
                            df = df.T[df.sum(axis=0) > 0].T
                            df = df[df.sum(axis=1) > 0]
                            df.to_hdf(dataFile + 'cut_%s.h5' % species, key='df', mode='a', complevel=4, complib='zlib')
                        else:
                            df = pd.read_hdf(dataFile + 'cut_%s.h5' % species, key='df')

                        DCS.prepare(df.astype(float))

                        df_expr = DCS.df_expr.loc[DCS.df_expr.index.intersection(np.unique(gEC23))][df_projection.columns]
                        seCluster = seCluster.reindex(df_expr.columns)

                        df_expr.columns = pd.MultiIndex.from_arrays([df_expr.columns.get_level_values('batch'),
                                                                     df_expr.columns.get_level_values('cell'),
                                                                     seCluster.values],
                                                                    names=['batch', 'cell', 'cluster'])
                        print(df_expr)

                        df_expr = df_expr.T[df_expr.sum(axis=0) > 0].T
                        df_expr = df_expr[df_expr.sum(axis=1) > 0].fillna(0.)
                        print(df_expr)

                        df_expr.sum(axis=0).to_excel('s0.xlsx')

                        #DCS.makeHeatmapGeneExpressionPlot(df_expr, selGenesHUGO)
                        DCS.makeHeatmapGeneExpressionPlot(df_expr, genes=gEC23, figsize = (6,9), fontsize=8, normalize=True, subtract=False)
            
                    ### Gene individual plots
                    if False:
                        #DCS.geneListFileName = 'geneLists/fibro_endo_epi_v2_Human.xlsx'
                        #markers = DCS.readMarkerFile().index.values.tolist()
                        #print('Markers:', len(markers))
                        #selGenesHUGO = [DCS.getHugoName(gene) for gene in markers + ['FLT4', 'VEGFA', 'VEGFB', 'VEGFC', 'VEGFD'] + tip + stalk]
                        #selGenesHUGO = [DCS.getHugoName(gene) for gene in ['ROBO4', 'ACTA2', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'WNT1', 'WNT2', 'WNT2B', 'WNT3', 'WNT3A', 'WNT4', 'WNT5A', 'WNT5B', 'WNT6', 'WNT7A', 'WNT7B', 'WNT8A', 'WNT8B', 'WNT9A', 'WNT9B', 'WNT10A', 'WNT10B', 'WNT11', 'WNT16']]

                        hemostasis_genes = pd.Series(np.loadtxt('hemostasis_genes.txt', dtype=str)).replace(Mouse_to_Human_HUGO_conversion).values.tolist()
                        selGenesHUGO = np.unique([DCS.getHugoName(gene) for gene in hemostasis_genes])
                        print(selGenesHUGO)

                        DCS.prepare(pd.read_hdf(dataFile + 'cut_%s.h5' % species, key='df').astype(float))

                        for gene in selGenesHUGO:
                            try:
                                expr = DCS.df_expr.xs(key=gene).loc[df_projection.columns].replace(0., np.nan).values
                                print(gene)
                                DCS.makeProjectionPlot(df_projection.values, reduce(expr), legend=False, labels=False, colorbar=True, suffix='by %s' % gene, dpi=150)
                            except:
                                pass

                    ### Plots by additional metadata
                    if False:
                        df_anno = df_anno.reindex(df_projection.columns.get_level_values('batch'), axis=0)

                        for col in ['medianOfGenesPerCell']:
                            DCS.makeProjectionPlot(df_projection.values, reduce(df_anno[col].values), legend=False, labels=False, colorbar=True, suffix='by %s cut' % col)

                        for col in ['isPrimaryAdult', 'isTumor', 'isCellLine']:
                            DCS.makeProjectionPlot(df_projection.values, df_anno[col].values, legend=True, labels=True, suffix='by %s cut' % col)

                    ### Plots by PnaglaoDB metadata
                    if False:
                        df_anno = df_anno.reindex(df_projection.columns.get_level_values('batch'), axis=0)

                        for col in ['Tissue', 'Species', 'Tissue detailed', 'Sequencing instrument', 'scRNA-seq protocol']:
                            DCS.makeProjectionPlot(df_projection.values, df_anno[col].values, legend=True, labels=True, suffix='by %s cut' % col)

                        for col in ['Sequenced cells', 'QC-passed fraction']:
                            DCS.makeProjectionPlot(df_projection.values, reduce(df_anno[col].values), legend=False, labels=False, colorbar=True, suffix='by %s cut' % col)

                        DCS.makeProjectionPlot(df_projection.values, df_projection.columns.get_level_values('batch').values, legend=False, labels=False, suffix='by batch cut')

                    ### Sankey plots
                    if False:
                        seTissue = pd.Series(index=df_projection.columns, data=df_anno['Tissue'].reindex(df_projection.columns.get_level_values('batch')).fillna(0.).values)
                        seTissueDetailed = pd.Series(index=df_projection.columns, data=df_anno['Tissue detailed'].reindex(df_projection.columns.get_level_values('batch')).fillna(0.).values)
                        seBatch = pd.Series(index=df_projection.columns, data=df_projection.columns.get_level_values('batch').values)
                        DCS.makeSankeyDiagram(DCS.getCountsDataframe(seBatch, seTissue), attemptSavingHTML=True, nameAppend=' Sankey batch-tissue')
                        DCS.makeSankeyDiagram(DCS.getCountsDataframe(seBatch, seTissueDetailed), attemptSavingHTML=True, nameAppend=' Sankey batch-tissueDet')
                        DCS.makeSankeyDiagram(DCS.getCountsDataframe(seTissueDetailed, seTissue), attemptSavingHTML=True, nameAppend=' Sankey tissueDet-tissue')

                        seCluster = pd.read_hdf(DCS.fileHDFpath, key='df_clusters').reindex(df_projection.columns)['cluster']
                        DCS.makeSankeyDiagram(DCS.getCountsDataframe(seCluster, seTissue), attemptSavingHTML=True, nameAppend=' Sankey cluster-tissue')
                        DCS.makeSankeyDiagram(DCS.getCountsDataframe(seCluster, seBatch), attemptSavingHTML=True, nameAppend=' Sankey cluster-batch')

                        seAnno = pd.read_hdf('data/seAnnotaionPanglaoDBsubsetEC.h5', key='seAnno/%s' % species)
                        sePanglaoCelltype = pd.Series(index=df_projection.columns, data=seAnno.reindex(df_projection.columns).fillna('Missing').values)
                        DCS.makeSankeyDiagram(DCS.getCountsDataframe(seTissue, sePanglaoCelltype), attemptSavingHTML=True, nameAppend=' Sankey tissue-celltype')
                        DCS.makeSankeyDiagram(DCS.getCountsDataframe(seTissueDetailed, sePanglaoCelltype), attemptSavingHTML=True, nameAppend=' Sankey tissueDet-celltype')
            
                # Find top-500 DE genes for each cluster (gene has to be expressed in at least 20% of cells of at least one cluster)
                if False:
                    cutoff = 1.
                    df_X = pd.read_hdf(dataFile + 'cut_%s.h5' % species, key='df').reindex(df_projection.columns, axis=1)
                    df_X.columns = pd.read_hdf(DCS.fileHDFpath, key='df_clusters').reindex(df_projection.columns)['cluster']
                    df_F = df_X.groupby(level=0, sort=True, axis=1).mean()
                    df_F = df_X.replace(0., np.nan).groupby(axis=1, level=0).agg('count') / df_X.fillna(0.).groupby(axis=1, level=0).agg('count')
                    print(df_F)

                    df_X = df_X.loc[df_F.max(axis=1) >= 0.20].groupby(level=0, sort=True, axis=1).mean()
                    print(df_X)

                    df_Z = (df_X - np.mean(df_X.values, axis=1)[:, None]) / np.std(df_X.values, axis=1)[:, None]
                    df_Z[df_Z < 0.3] = np.nan
                    print(df_Z)

                    sed = {col:df_Z[col].dropna().sort_values(ascending=False)[:500].index.values.tolist() for col in df_Z}

                    write(sed, 'Z dict 500 %s.gz' % species, jsonFormat=True)
                    print(sed)

                    print(df_Z.sum(axis=0))

                # Get Reactome sub-pathways hierarchy
                if True:
                    if True:
                        with open('rectomeGeneSets/ReactomePathways.gmt') as tempFile:
                            lines = [line[:-2].split('\t') for line in tempFile.readlines()]
                            lines = np.array([[line[1], line[0], line[2:]]  for line in lines])
                            reactomePathwayNames = dict(zip(lines.T[0], lines.T[1]))
                            write(reactomePathwayNames, 'reactome pathways names dict.gz', jsonFormat=True)

                        with open('rectomeGeneSets/ReactomePathwaysRelation.txt') as tempFile:
                            PathwaysRelation = np.array([line[:-2].split('\t') for line in tempFile.readlines()]).T
                            reactomePathwaysRelation = pd.Series(index=PathwaysRelation[0], data=PathwaysRelation[1])
                            write(reactomePathwaysRelation, 'reactome pathways relations dict')
                    else:
                        reactomePathwaysRelation = read('reactome pathways relations dict')
                        reactomePathwayNames = read('reactome pathways names dict.gz', jsonFormat=True)
            
                    reactomePathwayNames = pd.Series(reactomePathwayNames).apply(np.array)
                    print(reactomePathwayNames)
                    print(reactomePathwaysRelation)

                    #reactomePathwaysRelation = reactomePathwaysRelation[reactomePathwaysRelation.index.str.find('R-HSA-') >= 0]
                    ##reactomePathwaysRelation = pd.Series(index=reactomePathwaysRelation.values, data=reactomePathwaysRelation.index)
                    #reactomePathwaysRelation = reactomePathwaysRelation.groupby(level=0).agg('unique').apply(list).apply(lambda items: ', '.join(items))
                    #print(reactomePathwaysRelation)
                    #se_h = pd.Series(df.columns.get_level_values('pathwayID')).replace(reactomePathwaysRelation).values

                    pass

                # Combine fine cluster and make gene lists for Metascape analysis
                if False:
                    sed = read('Z dict 500 %s.gz' % species, jsonFormat=True)

                    df = pd.DataFrame(index=range(3*max(pd.Series(sed).apply(len))), columns=np.unique([key.split('.')[0] for key  in list(sed.keys())]))
                    for col in df.columns:
                        sel_genes = []
                        for sub in ['.0.0','.0.1','.0.2']:
                            try:
                                sel_genes.extend(sed[col + sub])
                            except:
                                pass
                        sel_genes = np.unique(sel_genes)

                        df[col] = pd.Series(sel_genes).reindex(df.index)
                    df.to_excel('forMetascape500.xlsx')

                # Analysis of cluster groups by their function
                if True:
                    # Reactome Analysis
                    if True:
                        if False:
                            sed = read('Z dict 500 %s.gz' % species, jsonFormat=True)

                            se = pd.Series(sed)
                            se = se[(se.apply(len) >= 10) & (se.apply(len) <= 1000)]
                            sed = se.to_dict()

                            from pyiomica.enrichmentAnalyses import ReactomeAnalysis, KEGGAnalysis, GOAnalysis, ExportEnrichmentReport, ExportReactomeEnrichmentReport

                            res = ReactomeAnalysis(sed)
                            res = {c:res[c].loc[res[c]['Entities FDR'] <= 0.05] for c in res.keys()}
                            res = {c:res[c] for c in res.keys() if len(res[c]) > 0}

                            write(res, 'temp_res %s' % species)
                            exit()

                            ExportReactomeEnrichmentReport(res, AppendString='Enrichment reports %s Reactome' % species)
                        else:

                            res = read('temp_res %s' % species)

                        # Export ordered summary file
                        if False:
                            df = pd.DataFrame(index=res.keys(), columns=np.unique(np.hstack([res[cluster].index.values for cluster in res]).flatten()))
                            df.columns = pd.MultiIndex.from_arrays([df.columns, pd.Series(df.columns).replace(reactomePathwayNames).values], names=['pathwayID', 'pathwayName'])

                            dictPathwaysInClusters = {cluster:res[cluster].index.values.tolist() for cluster in res}
                            for cluster in dictPathwaysInClusters:
                                df.loc[cluster, dictPathwaysInClusters[cluster]] = 1.

                            df.loc['sum'] = df.sum(axis=0)
                            df = df.sort_values(by='sum', ascending=False, axis=1).drop(['sum'], axis=0)
                            print(df)

                            def order(df):

                                data = df.fillna(0.).values

                                x = hierarchy.dendrogram(hierarchy.linkage(data, 'ward'), no_plot=True, get_leaves=True)['leaves']
                                y = hierarchy.dendrogram(hierarchy.linkage(data.T, 'ward'), no_plot=True, get_leaves=True)['leaves']

                                return df.iloc[x, y]

                            df = order(df).T
                            df.to_excel('v2 rep_500 %s.xlsx' % species)         
                            print(df)

                        if species == 'Mus musculus':
                            # We have identified 14 groups of clusters with dictinct functions, groups may partially overlap
                            groupClusters = ['9.0.0, 9.0.1, 9.0.2',
                                '13.0.0, 13.0.1, 13.0.2',
                                '6.0.1', 
                                '6.0.1, 6.0.2, 7.0.0, 7.0.2, 13.0.1', 
                                '2.0.2, 3.0.1, 5.0.0, 5.0.1, 6.0.0, 12.0.0, 12.0.1, 14.0.0', 
                                '5.0.0, 12.0.0, 12.0.1, 14.0.0', 
                                '1.0.0, 1.0.1, 1.0.2', 
                                '6.0.0, 14.0.0', 
                                '5.0.0, 14.0.0', 
                                '7.0.1', 
                                '2.0.0', 
                                '3.0.0, 3.0.1, 12.0.2', 
                                '0.0.0', 
                                '0.0.0, 0.0.2, 10.0.0, 10.0.1, 10.0.2']
                            groupFunction = ['Signalling',
                                'Phototransduction',
                                'Neuronal (muscle contraction)',
                                'Neuronal (signalling)',
                                'Extracellular matrix',
                                'Collagen',
                                'Metabolism',
                                'NOTCH3',
                                'NOTCH1',
                                'Cell cycle',
                                'BCR signalling',
                                'Surfactant metabolism',
                                'Costimulation',
                                'Immune']
                            dictGroupClustersToFunction = dict(zip(groupClusters, groupFunction))

                            if False:
                                dictGroupTotalCount = dict()
                                dictGroupPerTissueCount = dict()

                                seTissue = pd.Series(index=df_projection.columns, data=df_anno['Tissue'].reindex(df_projection.columns.get_level_values('batch')).fillna(0.).values)
                                seCluster = pd.read_hdf(DCS.fileHDFpath, key='df_clusters').reindex(df_projection.columns)['cluster']
                                dft = pd.concat([seTissue, seCluster], axis=1, sort=False)
                                dft.columns = ['tissue', 'cluster']
                                print(dft)

                                for group in groupClusters:
                                    cr = dict()
                                    for subc in group.split(', '):
                                        sel = dft['tissue'][dft['cluster'] == subc]
                                        u = np.unique(sel.values, return_counts=True)
                                        r = list(zip(u[0], u[1]))
                                        for tissue, count in r:
                                            if count >= 0:
                                                if tissue in cr.keys():
                                                    cr[tissue] += count
                                                else:
                                                    cr[tissue] = count

                                    cr = pd.Series(cr).sort_values(ascending=False)
                                    #print('\tSum:', cr.sum(), 'Cluster:', group, '\n', cr, '\n')

                                    dictGroupTotalCount.update({group:cr.sum()})
                                    dictGroupPerTissueCount.update({group:str(cr.to_dict()).replace("'", '').replace('}', '').replace('{', '')})

                            if False:
                                dfp = pd.read_excel('rep_500 Mus musculus.xlsx', sheet_name='cluster-sel-pathways', index_col=None, header=0)
                                dfp.columns = pd.MultiIndex.from_arrays([dfp.columns, pd.Series(dfp.columns).replace(dictGroupClustersToFunction).values], names=['groupClusters', 'groupFunction'])

                                dictHitsInPathways = {key:res[key]['Submitted entities found'].str.split(';').to_dict() for key in res}

                                cdf = pd.DataFrame(index=dfp.columns)

                                for group in dfp:
                                    clustersOfGroup = group[0].split(', ')
                                    pathwaysOfGroup = dfp[group].dropna().values.tolist()

                                    hitsOfGroup = []
                                    for cluster in clustersOfGroup:
                                        for pathway in pathwaysOfGroup:
                                            try:
                                                for gene in dictHitsInPathways[cluster][pathway]:
                                                    if not gene in hitsOfGroup:
                                                        hitsOfGroup.append(gene)
                                            except:
                                                pass

                                    cdf.loc[group, 'totalCount'] = dictGroupTotalCount[group[0]]

                                    for geneList in ['receptors', 'ligands', 'listTF', 'listUnlikelyTF', 'tip', 'stalk']:
                                        cdf.loc[group, geneList] = ','.join(sorted(list(set(hitsOfGroup).intersection(eval(geneList)))))

                                    cdf.loc[group, 'perTissueCount'] = dictGroupPerTissueCount[group[0]]
                                    cdf.loc[group, 'pathwaysOfGroup'] = ','.join(sorted(pathwaysOfGroup))

                                cdf = cdf.reset_index()
                                print(cdf)        

                                cdf.to_excel('cdf2.xlsx', index=False)

    # Correlation of receptors with all genes
    if False:
        # Combine samples into a single DataFrame
        if False:
            for species in ['Homo sapiens', 'Mus musculus']:
                DCS = DigitalCellSorter.DigitalCellSorter(dataName='PanglaoDB_EC_DCS', saveDir = 'DCS output/PanglaoDB_EC_DCS/%s/cut/'%species)
            
                cells = pd.read_hdf(DCS.fileHDFpath, key='df_projection').columns
                print('Cells:', len(cells))

                dataFile = 'DCS output/PanglaoDBendothelialAllv0_combined.h5'
                df = pd.read_hdf(dataFile, key='df').reindex(cells, axis=1)
                print('A', df)

                df = df[df.sum(axis=1) > 0]
                print('B', df)

                df = df.loc[:, df.sum(axis=0) > 0]
                print('C', df)

                df.to_hdf('DCS output/PanglaoDB_EC_byDCS_normed_and_filtered.h5', key=species, mode='a', complevel=4, complib='zlib')
                print('Done recording %s data to h5' % species)

        # Calculate fraction
        if False:
             for species in ['Homo sapiens', 'Mus musculus']:
                print('Reading normalized %s data from h5' % species, flush=True)
                df = pd.read_hdf('DCS output/PanglaoDB_EC_byDCS_normed_and_filtered.h5', key=species)
                print(df)

                df_fraction = df.replace(0., np.nan).groupby(axis=1, level='batch').agg('count') / df.fillna(0.).groupby(axis=1, level='batch').agg('count')
                print(df_fraction)
                df_fraction.to_hdf('DCS output/fraction_EC_gRec_by_species.h5', key=species, mode='a', complevel=4, complib='zlib')

        # Export expression of sel_genes
        if False:
             selGenes = np.unique(receptorsListHugo_2555 + gEC23)

             for species in ['Homo sapiens', 'Mus musculus']:
                print('Reading normalized %s data from h5' % species, flush=True)
                df = pd.read_hdf('DCS output/PanglaoDB_EC_byDCS_normed_and_filtered.h5', key=species)
                df = df.loc[df.index.intersection(selGenes)]
                print(df)

                df.to_hdf('DCS output/EC_expr_of_2555_and_EC23.h5', key=species, mode='a', complevel=4, complib='zlib')

        selGenes = np.unique(receptorsListHugo_2555 + gEC23)
        #corrFile = 'DCS output/correlation_EC_allRec_by_species_v9_2555.h5'
        #corrFile = 'DCS output/euclidean_EC_allRec_by_species_v10_2555.h5'
        corrFile = 'DCS output/euclidean_EC_allRec_by_species_v11_short_2555.h5'



        # Calculate euclidean distance of Human and Mouse separately
        if True:
            for species in ['Homo sapiens', 'Mus musculus']:
                print('Reading normalized %s %s data from h5' % ('EC', species), flush=True)
                df_corr = get_df_distance(pd.read_hdf('DCS output/PanglaoDB_EC_byDCS_normed_and_filtered_5k.h5', key=species), 
                                          invert=False, metric='euclidean', genes=selGenes, analyzeBy='batch', minSize=10)

                print('Recording to h5', flush=True)
                df_corr.to_hdf(corrFile, key=species, mode='a', complevel=4, complib='zlib')

        # Calculate corelation of Human and Mouse separately
        if False:
            for species in ['Homo sapiens', 'Mus musculus']:
                print('Reading normalized %s %s data from h5' % ('EC', species), flush=True)
                df_corr = get_df_distance(pd.read_hdf('DCS output/PanglaoDB_EC_byDCS_normed_and_filtered.h5', key=species), 
                                      genes=selGenes, analyzeBy='batch', minSize=10)

                print('Recording to h5', flush=True)
                df_corr.to_hdf(corrFile, key=species, mode='a', complevel=4, complib='zlib')

        # Test for number of optimal clusters
        if False:

            clusterAndExportCorrelations(corrFile, selGenes, 'EC', list(gEC23), TestNumClust=True)

        # Cluster correlations
        if False:

            clusterAndExportCorrelations(corrFile, selGenes, 'EC', list(gEC23), suff='v9')

        # Make dendrogram and heatmap
        if False:

            dendro(corrFile, selGenes, 'EC', list(gEC23), inhLoc = 18, n_clusters=10, saveDir='2555')

        # ROC calculation of various measures
        if False:

            rankPosNeg(corrFile, selGenes, 'EC', list(gEC23))


        if True:
            df = pd.read_excel('testing ROC AUC EC23 measures/compare Mouse Human top 50 metrics.xlsx',sheet_name='compare', index_col=None, header=[0,1])

            top = 50
            for top in [50]:
                print('Top: %s' % top)

                for metric in df.columns.levels[1]:
                    df_temp = df.xs(metric, axis=1, level=1).iloc[:top]
                    common = np.intersect1d(df_temp['mouse'].values, df_temp['human'].values)
                    common_h = np.intersect1d(df_temp['human'].values, gEC23)
                    common_m = np.intersect1d(df_temp['mouse'].values, gEC23)

                    print('Metric: %s,' % metric)
                    print('\tboth, any genes (%s):' % len(common), str(list(common)).replace("'", '').replace(']', '').replace('[', ''))
                    print('\tmouse (%s):' % len(common_m), str(list(common_m)).replace("'", '').replace(']', '').replace('[', ''))
                    print('\thuman (%s):' % len(common_h), str(list(common_h)).replace("'", '').replace(']', '').replace('[', ''))

                print()
