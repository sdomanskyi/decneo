import os
import numpy as np
import pandas as pd
from sergii_io import *
from string import Template

panglaoURL = Template('https://panglaodb.se/data_dl.php?sra=${SRA}&srs=${SRS}&filetype=R&datatype=readcounts')

def _selectCelltypeRecordH5(IndirName, fileHDF, samples, nMerge = 50):

    samples = np.unique(list(zip(df_sel.index.get_level_values(0).values, df_sel.index.get_level_values(1).values)), axis=0)

    for i, (SRA, SRS) in enumerate(samples):
        print('\n', i, 'of', len(samples), SRA, SRS)

        if KeyInStore('%s_%s' % (SRA, SRS), fileHDF):
            print('Key is already in h5 ifle')
            continue

        fileClusterMemberships = '%s%s.seurat_clusters.txt' % (SRA, '_' + SRS if SRS!='notused' else '')
        df_cell_to_cluster = pd.read_csv(os.path.join(IndirName, 'PanglaoDB', 'data', 'sample_clusters', fileClusterMemberships), delimiter=' ', index_col=0, header=None)[1]

        clustersToKeep = df_sel.xs(key=SRA, level=0).xs(key=SRS, level=0).index.values
        selectedCells = np.unique(df_cell_to_cluster[[(i in clustersToKeep) for i in df_cell_to_cluster.values]].index.values)

        dirName = os.path.join('..', '..', 'Downloads', 'var', 'www', 'html', 'SRA', 'SRA.final')

        df = readRDataFile('%s%s.sparse.RData' % (SRA, '_' + SRS if SRS!='notused' else ''), dirName)

        if df is None:
            print('Error reading df')

            continue

        df = df[selectedCells]
        print('Selected cells:', df.shape)

        print('Genes per cell:', (df!=0).sum(axis=0).median(), '\nReads per cell:', (df).sum(axis=0).median())
        nParts = int(np.ceil(len(df.columns) / nMerge))
        print(nParts)
        dfm = pd.concat([df[part].sum(axis=1) for part in np.array_split(df.columns, nParts)], axis=1, sort=False)
        dfm = dfm.loc[~dfm.index.duplicated(keep='first')]
        print('Genes per pseudo-cell:', (dfm!=0).sum(axis=0).median(), '\nReads per pseudo-cell:', (dfm).sum(axis=0).median())

        dfm.to_hdf(fileHDF, key='%s_%s' % (SRA, SRS), mode='a', complevel=4, complib='zlib')

    return

def getAnnotationsSummaryDf(dirName, saveToFile = True, printDf = False):

    df_metadata = pd.read_csv(os.path.join(dirName, 'PanglaoDB', 'data', 'metadata.txt'), index_col=[0, 1], header=None)
    df_metadata.index.names = ['SRA accession', 'SRS accession']
    df_metadata.columns = ['Tissue origin of the sample', 'scRNA-seq protocol', 'Species', 'Sequencing instrument', 'Number of expressed genes', 'Median number of expressed genes per cell', 'Number of cell clusters in this sample', 'Is the sample from a tumor? (1 true otherwise false)', 'Is the sample from primary adult tissue?', 'Is the sample from a cell line?']

    df_counts = pd.read_csv(os.path.join(dirName, 'PanglaoDB', 'data', 'counts.txt'), index_col=[0, 1], header=0).replace(0, np.nan)
    df_metadata = pd.concat([df_metadata, df_counts], axis=1, sort=False)
    df_metadata['Fraction of cells passed QC'] = df_metadata['Number of cells'] / df_metadata['Number of raw cells']

    df_metadata.sort_index(axis=1, inplace=True, ascending=False)

    df_cell_type_annotations = pd.read_csv(os.path.join(dirName, 'PanglaoDB', 'data', 'cell_type_annotations.txt'), index_col=[0, 1, 2], header=None)
    df_cell_type_annotations.index.names = ['SRA accession', 'SRS accession', 'Cluster index']
    df_cell_type_annotations.columns = ['Cell type annotation', 'P-value from Hypergeometric test', 'Adjusted p-value (BH)', 'Cell type Activity Score']

    df_clusters_to_number_of_cells = pd.read_csv(os.path.join(dirName, 'PanglaoDB', 'data', 'clusters_to_number_of_cells.txt'), index_col=[0, 1, 2], header=None)
    df_clusters_to_number_of_cells.index.names = ['SRA accession', 'SRS accession', 'Cluster index']
    df_clusters_to_number_of_cells.columns = ['Number of cells in cluster']
    df_cell_type_annotations = pd.concat([df_clusters_to_number_of_cells, df_cell_type_annotations], axis=1, sort=False)
    
    df_cell_type_annotations = df_cell_type_annotations.reindex(np.hstack([df_cell_type_annotations.columns, df_metadata.columns]), axis=1)
    df_cell_type_annotations.loc[:, df_metadata.columns] = df_metadata.loc[pd.MultiIndex.from_arrays([df_cell_type_annotations.index.get_level_values(0), df_cell_type_annotations.index.get_level_values(1)])].values

    df_cell_type_annotations = df_cell_type_annotations.replace('\\N', np.nan)

    for col in ['Cell type Activity Score',
                'Number of cell clusters in this sample',
                'Is the sample from primary adult tissue?',
                'Is the sample from a tumor? (1 true otherwise false)',
                'Is the sample from a cell line?']:
        df_cell_type_annotations[col] = df_cell_type_annotations[col].astype(float)

    celltypes = df_cell_type_annotations['Cell type annotation'].values
    df_cell_type_annotations = df_cell_type_annotations.loc[celltypes==celltypes]

    if saveToFile:
        if os.path.isfile(os.path.join(dirName, 'df_cell_type_annotations.xlsx')):
            print('Annotation summary file already exists')
        else:
            df_cell_type_annotations.to_excel(os.path.join(dirName, 'df_cell_type_annotations.xlsx'), merge_cells=False)

    if printDf:
        print(df_cell_type_annotations)

    return df_cell_type_annotations

def queryCounts(df_cell_type_annotations, MetadataDirName, RDataDirName):

    allSamples = np.unique(df_cell_type_annotations.droplevel(axis=0, level='Cluster index').index.values)

    result = [['SRA accession', 'SRS accession', 'Number of cells', 'Number of raw genes', 'Number of raw cells']]

    for i, (SRA, SRS) in enumerate(allSamples):
        print(i, 'of', len(allSamples), SRA, SRS)

        processedCells, rawGenes, rawCells = 0, 0, 0

        try:
            sourceFile = os.path.join(RDataDirName, '%s%s.sparse.RData.h5' % (SRA, '_' + SRS if SRS!='notused' else ''))
            rawGenes, rawCells = readRDataFile(sourceFile, returnSizeOnly=True)
        except Exception as exception:
            print(exception)

        try:
            processedCells = processedCells = df_cell_type_annotations.xs(level='SRA accession', key=SRA).xs(level='SRS accession', key=SRS)['Number of cells'].sum()
        except Exception as exception:
            print(exception)

        result.append([SRA, SRS, processedCells, rawGenes, rawCells])

    np.savetxt(os.path.join(MetadataDirName, 'counts.txt'), np.vstack(result), delimiter=',', fmt='%s')

    return

def extractPanglaoDBannotation(dirName, df_cell_type_annotations, SRA, SRS, saveFile = True):

    fileClusterMemberships = '%s%s.seurat_clusters.txt' % (SRA, '_' + SRS if SRS!='notused' else '')

    sourceFile = os.path.join(dirName, 'PanglaoDB', 'data', 'sample_clusters', fileClusterMemberships)
    df_cell_to_cluster = pd.read_csv(sourceFile, delimiter=' ', index_col=0, header=None)[1]
    df_cell_to_cluster.index.name = 'cell'

    panglaoClustersAnnotation = df_cell_type_annotations.xs(key=SRA, level=0).xs(key=SRS, level=0)['Cell type annotation'].to_dict()

    df_panglao_annotation = pd.concat([df_cell_to_cluster, df_cell_to_cluster.replace(panglaoClustersAnnotation)], axis=1, keys=['cluster', 'celltype'])

    if saveFile:
        df_panglao_annotation.to_hdf(os.path.join(dirName, 'PanglaoDBannotation', 'panglao_annotation_%s_%s.h5' % (SRA, SRS)), 
                                     key='df', mode='a', complevel=4, complib='zlib')

    return df_panglao_annotation

def getSRA(SRS, df_cell_type_annotations):

    return df_cell_type_annotations.xs(key=SRS, level='SRS accession').index.get_level_values('SRA accession')[0]

if __name__ == '__main__':

    dirName = os.path.join('dev', 'PanglaoDBdata', '')

    df_cell_type_annotations = getAnnotationsSummaryDf(dirName)

