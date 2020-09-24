import os
import numpy as np
import pandas as pd
from tools.io import read, write, KeyInStore
from rpy2.robjects import r
import urllib.request
from string import Template

def downloadFile(url, saveDir, saveName = None):

    if not os.path.exists(saveDir): 
        os.makedirs(saveDir)
    
    if saveName is None:
        saveName = url.strip('"').split('/')[-1:][0]

    path = os.path.join(saveDir, saveName)

    if os.path.isfile(path):
        print('File has been downloaded already')
    else:
        print('Downloading file:', url.strip('"'), end='\t', flush=True)

        try:
            urllib.request.urlretrieve(url.strip('"'), path)
            print('Done', flush=True)

        except Exception as exception:
            print(exception)

    return

def readPanglaoDBrDataFile(fileName, dirName, takeGeneSymbolOnly = True, maxSize = 2*10**9):

    fullPath = os.path.join(dirName, fileName)

    r['load'](fullPath)
    ls = np.array(r['ls']())
    #print('Variables in RData file:', ls)

    rSparseMatrix = r[ls[0]]

    dim = r('dim')(rSparseMatrix)
    
    if dim[0] * dim[1] > maxSize:
        print('Matrix size too large:', dim[0], dim[1])

        return None
    else:
        print('Matrix size:', dim[0], dim[1])

    rMatrix = r['as.matrix'](rSparseMatrix)

    df = pd.DataFrame(index=np.array(r['rownames'](rMatrix)),
                        columns=np.array(r['colnames'](rMatrix)),
                        data=np.array(rMatrix))

    if takeGeneSymbolOnly:
        df.index = df.index.str.split('_ENS', expand=True).get_level_values(0)

    return df

def getAnnotationsSummaryDf(dirName, saveToFile = True, printDf = True):

    df_metadata = pd.read_csv(os.path.join(dirName, 'PanglaoDB', 'data', 'metadata.txt'), index_col=[0, 1], header=None)
    df_metadata.index.names = ['SRA accession', 'SRS accession']
    df_metadata.columns = ['Tissue origin of the sample', 'scRNA-seq protocol', 'Species', 'Sequencing instrument', 'Number of expressed genes', 'Median number of expressed genes per cell', 'Number of cell clusters in this sample', 'Is the sample from a tumor? (1 true otherwise false)', 'Is the sample from primary adult tissue?', 'Is the sample from a cell line?']

    df_cell_type_annotations = pd.read_csv(os.path.join(dirName, 'PanglaoDB', 'data', 'cell_type_annotations.txt'), index_col=[0, 1, 2], header=None)
    df_cell_type_annotations.index.names = ['SRA accession', 'SRS accession', 'Cluster index']
    df_cell_type_annotations.columns = ['Cell type annotation', 'P-value from Hypergeometric test', 'Adjusted p-value (BH)', 'Cell type Activity Score']

    df_clusters_to_number_of_cells = pd.read_csv(os.path.join(dirName, 'PanglaoDB', 'data', 'clusters_to_number_of_cells.txt'), index_col=[0, 1, 2], header=None)
    df_clusters_to_number_of_cells.index.names = ['SRA accession', 'SRS accession', 'Cluster index']
    df_clusters_to_number_of_cells.columns = ['Number of cells']
    df_cell_type_annotations = pd.concat([df_clusters_to_number_of_cells, df_cell_type_annotations], axis=1, sort=False)
    
    df_cell_type_annotations = df_cell_type_annotations.reindex(np.hstack([df_cell_type_annotations.columns, df_metadata.columns]), axis=1)
    df_cell_type_annotations.loc[:, df_metadata.columns] = df_metadata.loc[pd.MultiIndex.from_arrays([df_cell_type_annotations.index.get_level_values(0), df_cell_type_annotations.index.get_level_values(1)])].values

    df_cell_type_annotations = df_cell_type_annotations.replace('\\N', np.nan)

    celltypes = df_cell_type_annotations['Cell type annotation'].values
    df_cell_type_annotations = df_cell_type_annotations.loc[celltypes==celltypes]

    if saveToFile:
        if os.path.isfile(os.path.join(dirName, 'df_cell_type_annotations.xlsx')):
            print('File already exists. Not overwriting it')
        else:
            df_cell_type_annotations.to_excel(os.path.join(dirName, 'df_cell_type_annotations.xlsx'))

    if printDf:
        print(df_cell_type_annotations)

    return df_cell_type_annotations

def selectAndRecordH5(IndirName, fileHDF, samples, nMerge = 50):

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

        #dirName = os.path.join('dev', 'PanglaoDBdata', '')
        #downloadFile('https://panglaodb.se/data_dl.php?sra=%s&srs=%s&filetype=R&datatype=readcounts' % (SRA, SRS), dirName, '%s_%s.sparse.RData' % (SRA, SRS))

        df = readPanglaoDBrDataFile('%s%s.sparse.RData' % (SRA, '_' + SRS if SRS!='notused' else ''), dirName)

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

def get_panglao_annotation_df(dirName, df_cell_type_annotations, SRA, SRS):

    fileClusterMemberships = '%s%s.seurat_clusters.txt' % (SRA, '_' + SRS if SRS!='notused' else '')
    df_cell_to_cluster = pd.read_csv(os.path.join(dirName, 'PanglaoDB', 'data', 'sample_clusters', fileClusterMemberships), delimiter=' ', index_col=0, header=None)[1]
    df_cell_to_cluster.index.name = 'cell'

    panglaoClustersAnnotation = df_cell_type_annotations.xs(key=SRA, level=0).xs(key=SRS, level=0)['Cell type annotation'].to_dict()

    df_panglao_annotation = pd.concat([df_cell_to_cluster, df_cell_to_cluster.replace(panglaoClustersAnnotation)], axis=1, keys=['cluster', 'celltype'])

    print(df_panglao_annotation)

    df_panglao_annotation.to_hdf(os.path.join(dirName, 'panglao_annotation_%s_%s.h5' % (SRA, SRS)), key='df', mode='a', complevel=4, complib='zlib')

    return df_panglao_annotation


if __name__ == '__main__':

    #dirName = os.path.join('')
    dirName = os.path.join('dev', 'PanglaoDBdata', '')
    panglaoURL = Template('https://panglaodb.se/data_dl.php?sra=${SRA}&srs=${SRS}&filetype=R&datatype=readcounts')

    if False:
        df_cell_type_annotations = getAnnotationsSummaryDf(dirName)

        #SRA, SRS = 'SRA713577', 'SRS3363004' # 3167 cells
        SRA, SRS = 'SRA550660', 'SRS2089637' # 9784 cells

        get_panglao_annotation_df(dirName, df_cell_type_annotations, SRA, SRS)
        downloadFile(panglaoURL.substitute(SRA=SRA, SRS=SRS), dirName, '%s_%s.sparse.RData' % (SRA, SRS))
        df = readPanglaoDBrDataFile('%s%s.sparse.RData' % (SRA, '_' + SRS if SRS!='notused' else ''), dirName)
        df.columns.names = ['cell']
        df.to_hdf(os.path.join(dirName, 'data_%s_%s.h5' % (SRA, SRS)), key='df', mode='a', complevel=4, complib='zlib')

        print(df)

    if True:
        import DigitalCellSorter

        #SRA, SRS = 'SRA713577', 'SRS3363004' # 3167 cells
        SRA, SRS = 'SRA550660', 'SRS2089637' # 9784 cells

        lists = ['CIBERSORT_LM22_7',
                 'CIBERSORT_LM22_14',
                 'CIBERSORT_LM22',
                 'cd_marker_handbook',
                 'PanglaoDB-Hs-44celltypes']

        geneListFileName = lists[4]

        DCS = DigitalCellSorter.DigitalCellSorter(dataName = 'PBMC',
                                                  geneListFileName = geneListFileName,
                                                  saveDir = os.path.join(os.path.dirname(__file__), 
                                                                         'DCS output', 
                                                                         'PBMC %s %s' % (SRA, SRS), 
                                                                         geneListFileName, ''))

        #DCS.prepare(pd.read_hdf(os.path.join(dirName, 'data_%s_%s.h5' % (SRA, SRS)), key='df'))
        #DCS.process()

        DCS.annotate()
        DCS.toggleMakeHistogramNullDistributionPlot = False
        DCS.visualize()

        df_panglao_annotation = pd.read_hdf(os.path.join(dirName, 'panglao_annotation_%s_%s.h5' % (SRA, SRS)), key='df')['celltype']
        df_panglao_annotation = df_panglao_annotation[df_panglao_annotation!=15]
        DCS.makeSankeyDiagram(DCS.getCountsDataframe(DCS.loadAnnotatedLabels(detailed=False), df_panglao_annotation)) 

    if False:
        df_cell_type_annotations = getAnnotationsSummaryDf(dirName)

        df_sel = df_cell_type_annotations.copy()
        df_sel = df_sel.loc[df_sel['Cell type annotation'].values == 'Dendritic cells']
        df_sel = df_sel.loc[df_sel['Species'].values == 'Homo sapiens']
        df_sel = df_sel.loc[df_sel['Is the sample from primary adult tissue?'].values.astype(str) == '1']
        df_sel = df_sel.loc[df_sel['Is the sample from a cell line?'].values.astype(str) == '0']
        df_sel = df_sel.loc[df_sel['Is the sample from a tumor? (1 true otherwise false)'].values.astype(str) == '0']
        df_sel = df_sel.loc[df_sel['Adjusted p-value (BH)'].values.astype(float) < 10**-2]

        selectAndRecordH5(dirName, 'dendritic.h5', df_sel)

        print('Reading h5')
        dfs = []
        with pd.HDFStore(os.path.join(dirName, 'dendritic.h5')) as hdf5file:
            keys = hdf5file.keys()

        for key in keys:
            print(key)
            df_temp = pd.read_hdf(os.path.join(dirName, 'dendritic.h5'), key=key)
            df_temp.columns = ['%s_%s' % (key.strip('/'), col) for col in list(range(df_temp.shape[1]))]
            dfs.append(df_temp)

        dfs = pd.concat(dfs, axis=1).fillna(0.)

        print('Before:', dfs.shape)
        dfs = dfs.loc[~dfs.index.duplicated(keep='first')]
        print('After:', dfs.shape)

        print(dfs)

        write(dfs, os.path.join(dirName, 'dfs_dendr'))

    if False:
        dirName = os.path.join('dev', 'PanglaoDBdata', '')

        #df_lnc = pd.read_csv('dev/gencode.v7.long_noncoding_RNAs.gtf', skiprows=[0,1,2,3,4], delimiter='\t')
        #df_lnc = df_lnc[df_lnc.columns[8]]
        #lnc = np.unique(df_lnc.str.split(';', expand=True)[4].str.split(' ', expand=True)[2].str.strip('"').str.split('.', expand=True)[0].values)
        #print(len(lnc))
        #temp = dfs.index.str.split('.', expand=True).get_level_values(0)
        #print(len(temp.intersection(lnc)))

        dfs = read(os.path.join(dirName, 'dfs_dendr'))
        print(dfs.shape)
        dfs = dfs.loc[~dfs.index.duplicated(keep='first')]
        print(dfs.shape)
        dfs.sort_index(inplace=True)

        #dfs.iloc[:, :1].to_excel('dev/x.xlsx'); exit()

        dfs = dfs.loc[(dfs>0).sum(axis=1)>0.25*dfs.shape[1]]
        print(dfs.shape)

        cutoff = int(0.25 * np.median(dfs.sum(axis=0)))
        print(cutoff)

        dfs = dfs.T.loc[(dfs.sum(axis=0)>cutoff).values].T

        s = dfs.sum(axis=0)

        lcutoff = int(0.5 * np.median(s))
        print('Lower:', lcutoff)

        ucutoff = int(2.5 * np.median(s))
        print('Upper:', ucutoff)

        dfs = dfs.T.loc[(s > lcutoff) & (s < ucutoff)].T

        import matplotlib.pyplot as plt
        dfs.sum(axis=0).hist(bins=50)
        plt.show()

        dfs.sort_index(inplace=True)


        #dfs.iloc[:, :3].to_excel('dev/x.xlsx')
        print(dfs.shape)




