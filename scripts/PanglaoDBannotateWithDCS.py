import scRegulation
from scRegulation.commonFunctions import *
from scRegulation import geneLists
import DigitalCellSorter

wfile = 'DCS output/PanglaoDBendothelialAllv0.h5'
cwfile = wfile[:-3] + '_combined.h5'
conffile = 'data/allConfidencePanglaoDB.h5'

if __name__ == '__main__':

    # Annotate all datasets of PanglaoDB
    if False:
        allFiles = np.loadtxt('PanglaoDBfilesSorted.txt', delimiter='\t', dtype=str)


        # For large datasets, where most cells are low quality according to Alona
        useQCcellsOnly = True
    
        confEC = []
        confEpC = []
        confFB = []
        for file in allFiles:
            if (len(file) > 16) and (file[-16:] == '.sparse.RData.h5'):
                print('\n\n\nProcessing', file, flush=True)

                batch = file[:-16]
                p = batch.split('_')
                SRA, SRS = p[0], 'notused' if len(p)==1 else p[1]

                saveDir = os.path.join('DCS output', 'PanglaoDB', batch, '')

                if os.path.isfile(os.path.join(saveDir, 'ColormapForCellTypes.txt')):
                    print('Directory already exists (%s). Not overriding it' % batch)
            
                    continue

                DCS = DigitalCellSorter.DigitalCellSorter(dataName=batch, precutQC=True, nClusters=12, doBatchCorrection=False, species='Human', mitochondrialGenesCutoffQC=1.75, saveDir=saveDir, geneListFileName = os.path.join('geneLists', 'fibro_endo_epi_v2_Human.xlsx'), updateConversionDictFile=False, layout='PCA')

                # Annotate cells of batch
                if False:
                    try:
                        df_expr = pd.read_hdf(os.path.join(RDataDirName, file), key='df')

                        if useQCcellsOnly:
                            df_expr = df_expr[pd.read_csv(os.path.join(MetadataDirName, 'PanglaoDB', 'data', 'sample_clusters', '%s%s.seurat_clusters.txt' % (SRA, '_' + SRS if SRS!='notused' else '')), delimiter=' ', index_col=0, header=None)[1].index.values]

                        df_expr.index = pd.Series(df_expr.index.values).replace(Mouse_to_Human_HUGO_conversion).values
                        df_expr = pd.concat([df_expr], keys=[batch], names=['batch'], axis=1)
                        df_expr.columns.names = ['batch', 'cell']
                        df_expr = df_expr.loc[~df_expr.index.duplicated(keep='first')]
                        df_expr = df_expr.T.loc[~df_expr.T.index.duplicated(keep='first')].T

                        df_expr0 = pd.read_hdf(os.path.join('data', 'Adult-Heart1_rmbatchdge.txt.gz.h5'), key='df')
                        df_expr0 = df_expr0.loc[~df_expr0.index.duplicated(keep='first')]
                        df_expr0 = df_expr0.T.loc[~df_expr0.T.index.duplicated(keep='first')].T
                
                        df_expr = pd.concat([df_expr, df_expr0], axis=1, sort=False)
                        df_expr = df_expr.T.loc[~df_expr.T.index.duplicated(keep='first')].T

                        DCS.excludedFromQC = df_expr.xs(key='AdultHeart_1', level='batch', axis=1, drop_level=False).columns

                        DCS.prepare(df_expr)
                        DCS.process()
                        DCS.annotate()

                        DCS.makeAnnotationResultsMatrixPlot()
                        DCS.makeMarkerExpressionPlot()
                        DCS.makeStackedBarplot()
                        DCS.makeProjectionPlotAnnotated()
                
                    except Exception as exception:
                        print('Error while annotating (%s)' % batch, exception)

                # Collect DCS endothelial cells into one hdf file
                if False:
                    try:
                        cells = DCS.getCells(celltype='Endothelial cells')

                        if not cells is None:
                            df_expr = DCS.getExprOfCells(cells)

                            df_expr = df_expr.xs(key=batch, axis=1, level='batch', drop_level=False)
                            df_expr = df_expr.loc[df_expr.sum(axis=1) > 0.]

                            columns = pd.MultiIndex.from_arrays([df_expr.columns.get_level_values('batch'), df_expr.columns.get_level_values('cell')])

                            df_expr = pd.DataFrame(data=df_expr.values, index=df_expr.index, columns=columns)
                            print('\tRemoved mixed-in cells and all-zero genes:', df_expr.shape)

                            print('\tSaving to hdf:', SRA, SRS)
                            df_expr.to_hdf(wfile, key=batch, mode='a', complevel=4, complib='zlib')

                    except Exception as exception:
                        print('Error while collecting Endothelial cells (%s)' % batch, exception)

                # Collect confidence EC
                if False:
                    try:
                        conf = pd.read_excel(os.path.join(DCS.saveDir, '%s_annotation.xlsx' % batch), sheet_name='z-scores', index_col='cluster', header=0)
                        se = pd.read_hdf(DCS.fileHDFpath, key='df_clusters')['cluster'].xs(key=batch, level='batch', drop_level=False)
                        seTemp = se.replace(conf['Endothelial cells'])
                        confEC.append(seTemp[seTemp > 0.])

                    except Exception as exception:
                        print('Error while annotating (%s)' % batch, exception)

                # Collect confidence EpC
                if False:
                    try:
                        conf = pd.read_excel(os.path.join(DCS.saveDir, '%s_annotation.xlsx' % batch), sheet_name='z-scores', index_col='cluster', header=0)
                        se = pd.read_hdf(DCS.fileHDFpath, key='df_clusters')['cluster'].xs(key=batch, level='batch', drop_level=False)
                        seTemp = se.replace(conf['Epithelial cells'])
                        confEpC.append(seTemp[seTemp > 0.])

                    except Exception as exception:
                        print('Error while annotating (%s)' % batch, exception)

                # Collect confidence FB
                if False:
                    try:
                        conf = pd.read_excel(os.path.join(DCS.saveDir, '%s_annotation.xlsx' % batch), sheet_name='z-scores', index_col='cluster', header=0)
                        se = pd.read_hdf(DCS.fileHDFpath, key='df_clusters')['cluster'].xs(key=batch, level='batch', drop_level=False)
                        seTemp = se.replace(conf['Fibroblasts'])
                        confFB.append(seTemp[seTemp > 0.])

                    except Exception as exception:
                        print('Error while annotating (%s)' % batch, exception)

        if len(confEC) > 0:
            confEC = pd.concat(confEC, axis=0, sort=False)
            confEC.to_hdf(conffile, key='Endothelial')
            print(confEC)

            confEpC = pd.concat(confEpC, axis=0, sort=False)
            confEpC.to_hdf(conffile, key='Epithelial')
            print(confEpC)

            confFB = pd.concat(confFB, axis=0, sort=False)
            confFB.to_hdf(conffile, key='Fibroblasts')
            print(confFB)
    
    # Combine samples to one DataFrame
    if False:
        keys = KeysOfStore(wfile)
        print('Reading %s keys of %s' % (len(keys), wfile), flush=True)

        dfs = []
        counts = 0
        for ikey, key in enumerate(keys):
            key = key[1:]

            print('\tReading:', ikey, 'of', len(keys), key, 'Count:', counts, flush=True)
            df_temp = pd.read_hdf(wfile, key=key).astype('float32')

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

        dfs.to_hdf(cwfile, key='df', mode='a', complevel=4, complib='zlib')

    # Meta-analysis
    if True:
        for species in ['Homo sapiens', 'Mus musculus'][:]:
            DCS = DigitalCellSorter.DigitalCellSorter(dataName='PanglaoDB_EC_DCS', verbose=2)
        
            # Process DCS-annotated EC
            if False:
                DCS.saveDir = 'results/DCS output/PanglaoDB_EC_DCS/%s/'%species

                # Un-cut EC clustering and analysis
                if True:
                    df = pd.read_hdf(cwfile, key='df')
                    df = df[df.columns[np.isin(df.columns.get_level_values('batch').values, df_anno[df_anno['Species'] == species].index)]]
                    df = df.T[df.sum(axis=0) > 0].T
                    df = df[df.sum(axis=1) > 0]
                    print(df)

                    DCS.prepare(df.astype(float))

                    del df

                    DCS.project()
                    DCS.nClusters = 15
                    DCS.cluster()

                # Select best endothelial cells
                if True:
                    df_projection = pd.read_hdf(DCS.fileHDFpath, key='df_projection')
                    valuesEndo = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Endothelial').reindex(df_projection.columns)
                    valuesEpi = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Epithelial').reindex(df_projection.columns)
                    valuesFibro = pd.read_hdf('data/allConfidencePanglaoDB.h5', key='Fibroblasts').reindex(df_projection.columns)
                    cutoff = 0.75
                    vEndo = ((valuesEndo.fillna(0.).values >= cutoff) & (valuesEpi.fillna(0.).values == 0.) & (valuesFibro.fillna(0.).values == 0.)).astype(bool)
                    print(len(vEndo), vEndo.sum(), vEndo.sum()/len(vEndo))
                    df_projection = df_projection[df_projection.columns[vEndo]]
                    print(df_projection.shape)

                DCS.saveDir = 'results/DCS output/PanglaoDB_EC_DCS/%s/cut/'%species

                # Cut EC clustering and analysis
                if True:
                    df = pd.read_hdf(cwfile, key='df')[df_projection.columns]
                    df = df.T[df.sum(axis=0) > 0].T
                    df = df[df.sum(axis=1) > 0]
                    print(df)

                    DCS.prepare(df.astype(float))

                    del df

                    DCS.project()
                    DCS.nClusters = 15
                    DCS.cluster()
                    DCS.makeProjectionPlotByClusters(suffix='by Clusters')
  
            # Make plots of filtered EC
            if True:
                DCS.saveDir='results/DCS output/PanglaoDB_EC_DCS/%s/cut/'%species

                df_projection = pd.read_hdf(DCS.fileHDFpath, key='df_projection')

                # Get PnaglaoDB metadata (df_anno)
                if True:
                    fields = ['Tissue origin of the sample', 'Species', 'Number of raw cells', 'Fraction of cells passed QC', 'Tissue origin of the sample', 'Sequencing instrument', 'scRNA-seq protocol', 'Median number of expressed genes per cell', 'Is the sample from primary adult tissue?', 'Is the sample from a tumor? (1 true otherwise false)', 'Is the sample from a cell line?']
                    df_anno = getPanglaoDBAnnotationsSummaryDf(MetadataDirName).droplevel('Cluster index', axis=0)[fields]
                    df_anno.columns = ['Tissue', 'Species', 'Sequenced cells', 'QC-passed fraction', 'Tissue detailed', 'Sequencing instrument', 'scRNA-seq protocol', 'medianOfGenesPerCell', 'isPrimaryAdult', 'isTumor', 'isCellLine']
                    df_anno = df_anno.loc[~df_anno.index.duplicated(keep='first')]
                    df_anno.index = ['%s%s' % (item[0], '_' + item[1] if item[1]!='notused' else '') for item in df_anno.index]

                    tissueDictPath = os.path.join(os.path.dirname(scRegulation.__file__), 'geneLists', 'PanglaoDBtissues.xlsx')
                    tissueDict = pd.read_excel(tissueDictPath, index_col='PanglaoDB tissue')['Starred'].str.replace('?', '')
                    tissueDict[tissueDict!=tissueDict] = tissueDict.index.values[tissueDict!=tissueDict]
                    df_anno['Tissue'] = df_anno['Tissue'].replace(tissueDict)

                # Sankey plots
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

                # Plots by PnaglaoDB metadata
                if True:
                    df_anno = df_anno.reindex(df_projection.columns.get_level_values('batch'), axis=0)

                    print('SRA:', len(np.unique(df_anno.index.str.split('_', expand=True).get_level_values(0))), 'SRS:', len(np.unique(df_anno.index.str.split('_', expand=True).get_level_values(1))))

                    makeBarplot(df_anno['Tissue'].values, DCS.saveDir, 'barplotByTissues')
                    DCS.makeProjectionPlot(df_projection.values, df_anno['Tissue'].values, legend=True, labels=True, suffix='by Tissue cut', rightShift=0.45) # 0.15 0.45

                    # Other quantities
                    if False:
                        for col in ['Tissue', 'isPrimaryAdult', 'isTumor', 'isCellLine', 'Species', 'Tissue detailed', 'Sequencing instrument', 'scRNA-seq protocol']:

                            DCS.makeProjectionPlot(df_projection.values, df_anno[col].values, legend=True, labels=True, suffix='by %s cut' % col, rightShift=0.45) # 0.15 0.45

                        for col in ['Sequenced cells', 'QC-passed fraction', 'medianOfGenesPerCell']:
                            DCS.makeProjectionPlot(df_projection.values, reduce(df_anno[col].values), legend=False, labels=False, colorbar=True, suffix='by %s cut' % col)

                        DCS.makeProjectionPlot(df_projection.values, df_projection.columns.get_level_values('batch').values, legend=False, labels=False, suffix='by batch cut')