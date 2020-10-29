import os
import numpy as np
import pandas as pd
from PanglaoDBannotation import *
import DigitalCellSorter
import h5py

import platform
if platform.system() == "Windows":
    RDataDirName = os.path.join('dev', 'PanglaoDBdata', '')
    MetadataDirName = os.path.join('dev', 'PanglaoDBdata', '')
else:
    RDataDirName = os.path.join('..', '..', 'Downloads', 'var', 'www', 'html', 'SRA', 'SRA.final')
    MetadataDirName = os.path.join('')

def readAllRData(RDataDirName, removeRDataFiles = True):

    filesToProcess = [file for file in os.listdir(RDataDirName) if file[-13:] == '.sparse.RData']
    
    for file in filesToProcess:
        fullFile = os.path.join(RDataDirName, file)

        try:
            if not os.path.isfile(fullFile + '.h5'):
                readRDataFile(fullFile + '.h5')

            if removeRDataFiles:
                print('Removing file', fullFile, flush=True)
                os.remove(fullFile)

        except Exception as exception:
            print(fullFile, '\n', exception, flush=True)

    return

def getCountsFromH5(RDataDirName):

    items = []
    for file in os.listdir(RDataDirName):
        if file[-9:] == '.RData.h5':
            with h5py.File(os.path.join(RDataDirName, file), 'r') as hdf5file:
                s = hdf5file['df/block0_values'].shape
                id = file.split('.sparse.RData.h5')[0].split('_')
                st = [id[0], 'notused' if len(id)==1 else id[1], s[0], s[1]]
                print(st)
                items.append(st)
    pd.DataFrame(items).set_index([0,1]).to_excel('h5_counts.xlsx', merge_cells = False)

    return

Mouse_to_Human_HUGO_conversion = read('geneLists/Mouse_to_Human_HUGO.gz', jsonFormat=True)

def tempFunc():
    
    DCSh = DigitalCellSorter.DigitalCellSorter(species='Human')
    DCSm = DigitalCellSorter.DigitalCellSorter(species='Mouse')

    conversion = read('Mouse_to_Human.gz', jsonFormat=True)
    conversion = pd.Series(conversion)
    conversion = pd.Series(index=DCSm.gnc.Convert(conversion.index.values.tolist(), sourceType='alias', targetType='hugo', returnUnknownString=False), data=DCSh.gnc.Convert(conversion.values.tolist(), sourceType='alias', targetType='hugo', returnUnknownString=False))
    conversion = conversion.loc[~conversion.index.duplicated(keep='first')].to_dict()

    write(conversion, 'Mouse_to_Human_HUGO.gz', jsonFormat=True)

    return

if __name__ == '__main__':

    df_cell_type_annotations = getAnnotationsSummaryDf(MetadataDirName)

    allFiles = np.loadtxt('PanglaoDBfilesSorted.txt', delimiter='\t', dtype=str)

    useQCcellsOnly = True

    # Collect confidence
    if False:
        allConfidenceEndo = []
        allConfidenceEpi = []
        allConfidenceFibro = []
        for ifile, file in enumerate(allFiles[:]):
            if (len(file) > 16) and (file[-16:] == '.sparse.RData.h5'):
                print('\n\n\nProcessing', ifile, file, flush=True)

                batch = file[:-16]
                p = batch.split('_')
                SRA, SRS = p[0], 'notused' if len(p)==1 else p[1]

                saveDir = os.path.join('DCS output', 'PanglaoDB', batch, '')

                DCS = DigitalCellSorter.DigitalCellSorter(dataName=batch, precutQC=True, nClusters=12, doBatchCorrection=False, species='Human', mitochondrialGenesCutoffQC=1.75, saveDir=saveDir, geneListFileName = os.path.join('geneLists', 'fibro_endo_epi_v2_Human.xlsx'), updateConversionDictFile=False, layout='PCA')

                try:
                    conf = pd.read_excel(os.path.join(DCS.saveDir, '%s_annotation.xlsx' % batch), sheet_name='z-scores', index_col='cluster', header=0)
                    se = pd.read_hdf(DCS.fileHDFpath, key='df_clusters')['cluster'].xs(key=batch, level='batch', drop_level=False)

                    seEndo = se.replace(conf['Endothelial cells'])
                    seEpi = se.replace(conf['Epithelial cells'])
                    seFibro = se.replace(conf['Fibroblasts'])

                    allConfidenceEndo.append(seEndo[seEndo > 0.])
                    allConfidenceEpi.append(seEpi[seEpi > 0.])
                    allConfidenceFibro.append(seFibro[seFibro > 0.])

                except Exception as exception:
                    print('Error while annotating (%s)' % batch, exception)

        allConfidenceEndo = pd.concat(allConfidenceEndo, axis=0, sort=False)
        print(allConfidenceEndo)
        allConfidenceEndo.to_hdf('data/allConfidencePanglaoDB.h5', key='Endothelial')

        allConfidenceEpi = pd.concat(allConfidenceEpi, axis=0, sort=False)
        print(allConfidenceEpi)
        allConfidenceEpi.to_hdf('data/allConfidencePanglaoDB.h5', key='Epithelial')

        allConfidenceFibro = pd.concat(allConfidenceFibro, axis=0, sort=False)
        print(allConfidenceFibro)
        allConfidenceFibro.to_hdf('data/allConfidencePanglaoDB.h5', key='Fibroblasts')

    # Annotate and collect cells
    if False:
        batch1 = allFiles[0:200]
        batch2 = allFiles[200:400]
        batch3 = allFiles[400:600]
        batch4 = allFiles[600:800]
        batch5 = allFiles[800:1000]
        batch6 = allFiles[1000:1100]
        batch7 = allFiles[1100:1200]
        batch8 = allFiles[1200:1300]
        batch9 = allFiles[1300:]

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

                if False:
                    try:
                        DCS.makeSankeyDiagram(DCS.getCountsDataframe(DCS.loadAnnotatedLabels(detailed=True), DCS.loadAnnotatedLabels(infoType='batch')), nameAppend='_Sankey_diagram_batch')

                        PanglaoDBannotation = extractPanglaoDBannotation(MetadataDirName, df_cell_type_annotations, SRA, SRS)
                        df_counts = DCS.getCountsDataframe(DCS.loadAnnotatedLabels(), PanglaoDBannotation)
                        print(df_counts)
                        DCS.makeSankeyDiagram(df_counts, attemptSavingHTML=True, nameAppend='_Sankey_diagram_PanglaoAnnotation') 

                    except Exception as exception:
                        print('Error while comparing to PanglaoDB alona annotation (%s)' % batch, exception)

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
                            df_expr.to_hdf(os.path.join('DCS output', 'PanglaoDBendothelialAllv0.h5'), key=batch, mode='a', complevel=4, complib='zlib')

                    except Exception as exception:
                        print('Error while collecting Endothelial cells (%s)' % batch, exception)

                if False:
                    try:
                        cells = DCS.getCells(celltype='Epithelial cells')

                        if not cells is None:
                            df_expr = DCS.getExprOfCells(cells)

                            df_expr = df_expr.xs(key=batch, axis=1, level='batch', drop_level=False)
                            df_expr = df_expr.loc[df_expr.sum(axis=1) > 0.]

                            columns = pd.MultiIndex.from_arrays([df_expr.columns.get_level_values('batch'), df_expr.columns.get_level_values('cell')])

                            df_expr = pd.DataFrame(data=df_expr.values, index=df_expr.index, columns=columns)
                            print('\tRemoved mixed-in cells and all-zero genes:', df_expr.shape)

                            print('\tSaving to hdf:', SRA, SRS)
                            df_expr.to_hdf(os.path.join('DCS output', 'PanglaoDBepithelialAllv0.h5'), key=batch, mode='a', complevel=4, complib='zlib')

                    except Exception as exception:
                        print('Error while collecting Epithelial cells (%s)' % batch, exception)

                if False:
                    try:
                        cells = DCS.getCells(celltype='Fibroblasts')

                        if not cells is None:
                            df_expr = DCS.getExprOfCells(cells)

                            df_expr = df_expr.xs(key=batch, axis=1, level='batch', drop_level=False)
                            df_expr = df_expr.loc[df_expr.sum(axis=1) > 0.]

                            columns = pd.MultiIndex.from_arrays([df_expr.columns.get_level_values('batch'), df_expr.columns.get_level_values('cell')])

                            df_expr = pd.DataFrame(data=df_expr.values, index=df_expr.index, columns=columns)
                            print('\tRemoved mixed-in cells and all-zero genes:', df_expr.shape)

                            print('\tSaving to hdf:', SRA, SRS)
                            df_expr.to_hdf(os.path.join('DCS output', 'PanglaoDBfibroblastsAllv0.h5'), key=batch, mode='a', complevel=4, complib='zlib')

                    except Exception as exception:
                        print('Error while collecting Fibroblasts (%s)' % batch, exception)
