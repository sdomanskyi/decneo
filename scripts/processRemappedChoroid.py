from scRegulation.commonFunctions import *
import DigitalCellSorter

def read10xSaveDonorChoroid(donor, inDir = ''):

    dir = inDir + donor + '/outs/filtered_feature_bc_matrix/'

    df_data = pd.DataFrame(data=scipy.io.mmread(dir + 'matrix.mtx.gz').tocsr().todense(), 
                           index=pd.read_csv(dir + 'features.tsv.gz', compression='gzip', sep='\t', index_col=1, header=None).index, 
                           columns=pd.read_csv(dir + 'barcodes.tsv.gz', compression='gzip', sep='\t', index_col=0, header=None).index)

    df_data = df_data.loc[df_data.sum(axis=1) > 0.]
    df_data = df_data.T.loc[df_data.T.sum(axis=1) > 0.].T

    print(df_data)

    df_data.to_hdf('Voigt_choroid_4_5_6_7_remapped_bySD_counts.h5', key=donor, mode='a', complevel=4, complib='zlib')

    return

def collectRemappedChoroid():

    if False:
        for i in [4, 5, 6, 7]:
            for location in ['macula', 'peripheral']:
                read10xSaveDonorChoroid('%s_donor_%s' % (location, i))

        dfs = []
        for i in [4, 5, 6, 7]:
            for location in ['macula', 'peripheral']:
                donor = '%s_donor_%s' % (location, i)
                df_temp = pd.read_hdf('Voigt_choroid_4_5_6_7_remapped_bySD_counts.h5', key=donor)
                df_temp.columns.names = ['cell']
                df_temp.index.names = ['gene']
                df_temp = pd.concat([df_temp], axis=1, keys=[donor], names=['batch']).reorder_levels(['batch', 'cell'], axis=1)
                df_temp = df_temp.loc[~df_temp.index.duplicated(keep='first')]
                df_temp = df_temp.T.loc[~df_temp.columns.get_level_values('cell').duplicated(keep='first')].T
                print(df_temp)
                dfs.append(df_temp)

        dfs = pd.concat(dfs, axis=1).fillna(0).astype(int)
        print(dfs)
        dfs.to_hdf('Voigt_choroid_4567_remapped_bySD_counts.h5', key='df', mode='a', complevel=4, complib='zlib')

        os.remove('Voigt_choroid_4_5_6_7_remapped_bySD_counts.h5')

    if False:
        dfs = []
        dir = '../data/byAV/'
        for file in os.listdir(dir):
            df_temp = pd.read_csv(dir + file, compression='gzip', sep='\t', header=0, index_col=[0,1,2,3,4,5]).T.droplevel(['final_cluster_labels', 'donor', 'region', 'age'], axis=1)
            df_temp.columns.names= ['cell', 'celltype']
            df_temp = pd.concat([df_temp], axis=1, keys=[file.strip('_CD31_pos_COUNTS.tsv.gz')], names=['batch']).reorder_levels(['batch', 'cell', 'celltype'], axis=1)
            df_temp = df_temp.loc[~df_temp.index.duplicated(keep='first')]
            df_temp = df_temp.T.loc[~df_temp.columns.get_level_values('cell').duplicated(keep='first')].T
            print(df_temp)
            dfs.append(df_temp)

        dfs = pd.concat(dfs, axis=1).fillna(0).astype(int)
        print(dfs)
        dfs.to_hdf('Voigt_choroid_4567_byAV_counts.h5', key='df', mode='a', complevel=4, complib='zlib')

    return

def makeGeneProjectionPlot(genes, dataName='dataName', saveDir=''):

    DCS = DigitalCellSorter.DigitalCellSorter(dataName=dataName, saveDir=saveDir, layout='TSNE')
    df_projection = pd.read_hdf(DCS.fileHDFpath, key='df_projection')

    DCS.loadExpressionData()

    if type(genes) is str:
        gene = [genes]

    for gene in genes:
        try:
            values = reduce(DCS.df_expr.loc[gene].values.astype(float))
            values[values == 0.] = np.nan
            DCS.makeProjectionPlot(df_projection.values, values, legend=False, labels=False, colorbar=True, suffix='by_%s' % gene, rightShift=0.2)
        except:
            print('ERROR:', gene)
    return

if __name__ == '__main__':

    if False:
        DCS = DigitalCellSorter.DigitalCellSorter(saveDir='remappedChoroid4567_DCS/')

        df_SD = pd.read_hdf('Voigt_choroid_4567_remapped_bySD_counts.h5', key='df')
        df_SD.index = DigitalCellSorter.DigitalCellSorter().gnc.Convert(df_SD.index.values.tolist(), 'alias', 'hugo', returnUnknownString=False)
        df_SD = df_SD[df_SD.columns[(((df_SD>0).sum(axis=0)>=1101) & (df_SD.sum(axis=0)>=1775))]]
        donorsLevel = df_SD.columns.levels[0].copy()
        #df_SD = df_SD.iloc[0].reset_index().set_index('cell')
        #df_SD.index = df_SD['batch'].index.str.split('-', expand=True).get_level_values(0).str.split('_', expand=True).get_level_values(0)
        df_SD = df_SD.loc[~df_SD.index.duplicated(keep=False)]
        print(df_SD)

        df_AV = pd.read_hdf('Voigt_choroid_4567_byAV_counts.h5', key='df').droplevel('celltype', axis=1)
        df_AV.columns = df_AV.columns.set_levels(donorsLevel, level=0)
        df_AV.index = DigitalCellSorter.DigitalCellSorter().gnc.Convert(df_AV.index.values.tolist(), 'alias', 'hugo', returnUnknownString=False)
        #df_AV = df_AV.iloc[0].reset_index().set_index('cell')
        #df_AV.index = df_AV['batch'].index.str.split('-', expand=True).get_level_values(0).str.split('_', expand=True).get_level_values(0)
        df_AV = df_AV.loc[~df_AV.index.duplicated(keep=False)]
        print(df_AV)

        se_SD = (df_SD>0).mean(axis=1)
        se_SD = se_SD[se_SD>0]

        se_AV = (df_AV>0).mean(axis=1)
        se_AV = se_AV[se_AV>0]

        print(se_SD)
        print(se_AV)

        se_AV[se_AV.index.difference(se_SD.index)].sort_values(ascending=False).to_excel('AV_only_H.xlsx')
        se_SD[se_SD.index.difference(se_AV.index)].sort_values(ascending=False).to_excel('SD_only_H.xlsx')
        pd.concat([se_SD[se_SD.index.intersection(se_AV.index)], se_AV[se_SD.index.intersection(se_AV.index)]], axis=1, sort=False).sort_values(by=0, ascending=False).to_excel('SD_and_AV_H.xlsx')

        DCS.makeSankeyDiagram(DCS.getCountsDataframe(df_SD['batch'], df_AV['batch']))

        exit()

    DCS = DigitalCellSorter.DigitalCellSorter(saveDir='results/remappedChoroid4567_DCS/',
                                              precutQC=True, nClusters=10, useUnderlyingNetwork=True, doFineClustering=True,
                                              minSubclusterSize=25, minimumNumberOfMarkersPerCelltype=3, thresholdForUnknown=0.2,
                                              splitFineClusters=True, subSplitSize=50, doBatchCorrection=False, zScoreCutoff=0.2,
                                              geneListFileName='geneLists/eye_correlated_markers.xlsx', verbose=1)

    if False:
        DCS.prepare(pd.read_hdf('Voigt_choroid_4567_remapped_bySD_counts.h5', key='df'))
        DCS.process()
        DCS.annotate()

        DCS.makeProjectionPlotByBatches()
        DCS.makeStackedBarplot()
        DCS.makeHopfieldLandscapePlot(plot3D=True, attemptSavingHTML=True)
        DCS.makeHopfieldLandscapePlot(plot3D=False, reuseData=True)
        DCS.makeHopfieldPCplot(trPath=DCS.saveDir + 'Hopfiled_mode_1/')

        DCS.makeProjectionPlotsQualityControl()
        DCS.makeMarkerExpressionPlot()

        DCS.makeProjectionPlotAnnotated(legend=False)

    if True:
        se = DCS.loadAnnotatedLabels()

        makeBarplot(se.values, DCS.saveDir, 'barplotBy_celltype')
        makeBarplot(se[se!='Failed QC'].index.get_level_values('batch'), DCS.saveDir, 'barplotBy_batches')
        makeBarplot(se[se == 'Endothelial'].index.get_level_values('batch'), DCS.saveDir, 'barplotEC_By_batches')


    if False:
        df_projection = pd.read_hdf(DCS.fileHDFpath, key='df_projection')
        DCS.makeProjectionPlot(df_projection.values, df_projection.columns.get_level_values('batch'), legend=True, labels=True, suffix='byBatches', rightShift = 0.4)
        se_celltype = pd.read_hdf(DCS.fileHDFpath, key='df_markers_expr').iloc[0].reset_index().set_index(['batch', 'cell'])['label'].str.split(' #', expand=True)[0]
        DCS.makeProjectionPlot(df_projection.values, se_celltype.reindex(df_projection.columns).values, legend=True, labels=True, suffix='byCelltype', rightShift = 0.4)

        makeGeneProjectionPlot(['NRP1', 'ITGA6', 'PECAM1', 'KDR', 'LIFR', 'IL6ST', 'ROBO4', 'PLXND1', 'OSMR'], saveDir='remappedChoroid4567_DCS/')

    # ?
    if False:
        se = DCS.loadAnnotatedLabels()
        print(se)
        se = se[se != 'Failed QC']
        se.to_hdf('Voigt_choroid_4567_remapped_bySD_annotation_labels_13700cells.h5', key='df', mode='a', complevel=4, complib='zlib')
        print(se)
        exit()

        DCS.loadExpressionData()
        DCS.df_expr.droplevel('cluster', axis=1).to_hdf('Voigt_choroid_4567_remapped_bySD_DCS_all.h5', key='df', mode='a', complevel=4, complib='zlib')
        exit()

        df_expr = DCS.getExprOfCells(DCS.getCells('Endothelial')).droplevel('cluster', axis=1)
        print(df_expr)
        df_expr.to_hdf('Voigt_choroid_4567_remapped_bySD_DCS_EC.h5', key='df', mode='a', complevel=4, complib='zlib')

        se = DCS.loadAnnotatedLabels()
        se.index = se.index.set_levels([cell.split('-')[0] for cell in se.index.levels[1]], level=1)
        print(se)

        df_AV = pd.read_hdf('Voigt_choroid_4567_byAV_counts.h5', key='df')
        df_AV.columns = df_AV.columns.set_levels(se.index.levels[0], level=0)
        df_AV.index = DigitalCellSorter.DigitalCellSorter().gnc.Convert(df_AV.index.values.tolist(), 'alias', 'hugo', returnUnknownString=False)
        df_AV = df_AV.loc[~df_AV.index.duplicated(keep=False)]

        df_AV.columns = pd.MultiIndex.from_arrays([df_AV.columns.get_level_values('batch'), df_AV.columns.get_level_values('cell').str.split('_', expand=True).get_level_values(0), df_AV.columns.get_level_values('celltype')], names=['batch', 'cell', 'celltype'])

        df_AV = df_AV.iloc[0].reset_index().set_index(['batch', 'cell'])['celltype']
        print(df_AV)

        DCS.makeSankeyDiagram(DCS.getCountsDataframe(se, df_AV), nameAppend='v2', attemptSavingHTML=True)
