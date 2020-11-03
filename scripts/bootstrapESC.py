from decneo.commonFunctions import *
from decneo.analysisPipeline import Analysis, process

def prepareInputData_human_McCracken():

    dir = '/mnt/research/piermarolab/data/3. McCracken GSE131736/pd/'

    def readOne(fileName, celltype, batch):

        df = pd.read_hdf(dir + fileName, key='df')
        df = pd.concat([df], keys=[celltype], axis=1, sort=False)
        df = pd.concat([df], keys=[batch], axis=1, sort=False)
        df.columns.names = ['batch', 'celltype', 'cell']
        df = df.reorder_levels(['batch', 'cell', 'celltype'], axis=1)

        return df

    df = pd.concat([readOne('other/GSM3814885_h9_day0.h5', 'day0_hESC', 'GSE131736'),
                    readOne('GSM3814888_h9_day8_rep1.h5', 'day8_mesodermal', 'GSE131736').sample(n=3000, axis=1),
                    readOne('GSM3814889_h9_day8_rep2.h5', 'day8_mesodermal', 'GSE131736').sample(n=3000, axis=1),
                    readOne('GSM3814890_h9_day8_rep3.h5', 'day8_mesodermal', 'GSE131736').sample(n=3000, axis=1)],
                    axis=1, sort=False).fillna(0.).astype(int)
    print(df)

    import DigitalCellSorter
    Convert = DigitalCellSorter.DigitalCellSorter(matplotlibMode='TkAgg').gnc.Convert
    df.index = Convert(df.index.values.tolist(), 'alias', 'hugo', returnUnknownString=False)
    df = df.loc[~df.index.duplicated(keep='first')]
    df = df.T.loc[~df.T.index.get_level_values('cell').duplicated(keep='first')].T
    df = df.astype(float)
    df = df[df.sum(axis=1) > 0].apply(lambda cell: cell * 10000. / np.sum(cell), axis=0)
    df = np.log2(df.replace(0., np.min(df.values[df.values > 0.])))
    df -= np.min(df.values)
    df = df[np.std(df, axis=1) / np.mean(np.std(df.values)) > 0.01]
        
    isEC = df.columns.get_level_values('celltype').values == 'day0_hESC'
    df = df.droplevel('celltype', axis=1).astype(float)
    df_EC = df[df.columns[isEC]]
    df_other = df[df.columns[~isEC]]

    nBatches = 10
    EC_batches = np.hstack([v + '_' + str(i) for i, v in enumerate(np.array_split(df_EC.columns.get_level_values('batch').values, nBatches))])
    other_batches = np.hstack([v + '_' + str(i) for i, v in enumerate(np.array_split(df_other.columns.get_level_values('batch').values, nBatches))])
    np.random.shuffle(EC_batches)
    np.random.shuffle(other_batches)

    df_EC.columns = pd.MultiIndex.from_arrays([EC_batches, df_EC.columns.get_level_values('cell')], names=['batch', 'cell'])
    df_other.columns = pd.MultiIndex.from_arrays([other_batches, df_other.columns.get_level_values('cell')], names=['batch', 'cell'])

    print(df_EC)
    print(df_other)

    return df_EC, df_other

def prepareInputData_mouse_Han():

    dir = '/mnt/research/piermarolab/data/5. MCA/MCA1.0/rmbatch_dge/'

    def readOne(fileName, celltype, batch):

        df = pd.read_hdf(dir + fileName, key='df').droplevel('batch', axis=1)
        df = pd.concat([df], keys=[celltype], axis=1, sort=False)
        df = pd.concat([df], keys=[batch], axis=1, sort=False)
        df.columns.names = ['batch', 'celltype', 'cell']
        df = df.reorder_levels(['batch', 'cell', 'celltype'], axis=1)

        return df

    df = pd.concat([readOne('EmbryonicStemCells_rm.batch_dge.txt.gz.h5', 'mESC', 'SRA638923'),
                    readOne('EmbryonicMesenchymeE14.5_rm.batch_dge.txt.gz.h5', 'other', 'SRA638923')],
                    axis=1, sort=False).fillna(0.).astype(int)

    import DigitalCellSorter
    Convert = DigitalCellSorter.DigitalCellSorter(matplotlibMode='TkAgg').gnc.Convert
    df.index = Convert(df.index.values.tolist(), 'alias', 'hugo', returnUnknownString=False)
    df = df.loc[~df.index.duplicated(keep='first')]
    df = df.T.loc[~df.T.index.get_level_values('cell').duplicated(keep='first')].T
    df = df.astype(float)
    df = df[df.sum(axis=1) > 0].apply(lambda cell: cell * 10000. / np.sum(cell), axis=0)
    df = np.log2(df.replace(0., np.min(df.values[df.values > 0.])))
    df -= np.min(df.values)
    df = df[np.std(df, axis=1) / np.mean(np.std(df.values)) > 0.01]
        
    isEC = df.columns.get_level_values('celltype').values == 'mESC'
    df = df.droplevel('celltype', axis=1).astype(float)
    df_EC = df[df.columns[isEC]]
    df_other = df[df.columns[~isEC]]

    nBatches = 10
    EC_batches = np.hstack([v + '_' + str(i) for i, v in enumerate(np.array_split(df_EC.columns.get_level_values('batch').values, nBatches))])
    other_batches = np.hstack([v + '_' + str(i) for i, v in enumerate(np.array_split(df_other.columns.get_level_values('batch').values, nBatches))])
    np.random.shuffle(EC_batches)
    np.random.shuffle(other_batches)

    df_EC.columns = pd.MultiIndex.from_arrays([EC_batches, df_EC.columns.get_level_values('cell')], names=['batch', 'cell'])
    df_other.columns = pd.MultiIndex.from_arrays([other_batches, df_other.columns.get_level_values('cell')], names=['batch', 'cell'])

    print(df_EC)
    print(df_other)

    return df_EC, df_other


if __name__ == '__main__':

    wdir = '/mnt/research/piermarolab/Sergii/results/ESC/'

    process(*prepareInputData_human_McCracken(), *prepareInputData_mouse_Han(),
            wdir + 'McCracken_hESC_vs_day8/', wdir + 'Han_mESC_vs_mesenchyme/', 
            nCPUs=8 if platform.system()=="Windows" else 10, parallelBootstrap=True,
            genesOfInterest=TF, knownRegulators=TFmarkers, exprCutoff1=0.05, exprCutoff2=0.005, perEachOtherCase=True)