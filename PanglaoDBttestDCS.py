import os
import sys
import numpy as np
import pandas as pd
import DigitalCellSorter
from commonFunctions import *

if __name__ == '__main__':

    if False:
        ECcells1 = pd.Series(index=pd.read_hdf('results/DCS output/PanglaoDB_EC_byDCS_normed_and_filtered.h5', key='Homo sapiens').columns, data='Homo sapiens')
        ECcells2 = pd.Series(index=pd.read_hdf('results/DCS output/PanglaoDB_EC_byDCS_normed_and_filtered.h5', key='Mus musculus').columns, data='Mus musculus')
        ECboth = pd.concat([ECcells1, ECcells2], axis=0, sort=False)
        ECboth.to_hdf('results/DCS output/PanglaoDB_ECcells.h5', key='df', mode='a', complevel=4, complib='zlib')
        print(ECboth)
    else:

        ECboth = pd.read_hdf('results/DCS output/PanglaoDB_ECcells.h5', key='df')

    if False:

        fileName = 'results/DCS output/ttest_PanglaoEC.h5'

        dir = 'results/DCS output/PanglaoDB/'
        for i, batch in enumerate(os.listdir(dir)[:]):
            print('Processing', i, batch, end='\t')
            DCS = DigitalCellSorter.DigitalCellSorter(dataName=batch, saveDir=dir + batch)

            try:
                ECcells = ECboth.xs(key=DCS.dataName, axis=0, level='batch', drop_level=False).index

                if len(ECcells) >= 10:
                    DCS.loadExpressionData()
                    df_expr = DCS.df_expr.xs(key=DCS.dataName, axis=1, level='batch', drop_level=False)
                    columns = pd.MultiIndex.from_arrays([df_expr.columns.get_level_values('batch'), df_expr.columns.get_level_values('cell')])
                    df_expr = pd.DataFrame(data=df_expr.values, index=df_expr.index, columns=columns)

                    df_EC = df_expr[ECcells]
                    df_other = df_expr[df_expr.columns.difference(df_EC.columns)]
                
                    df_ttest = pd.DataFrame(index=df_EC.index, columns=['statistic', 'pvalue'])
                    ttest = scipy.stats.ttest_ind(df_EC.values, df_other.values, axis=1)
                    df_ttest['statistic'] = ttest[0]
                    df_ttest['pvalue'] = ttest[1]
                    df_ttest = df_ttest.sort_values('statistic', ascending=False).dropna()
                    df_ttest.to_hdf(fileName, key=DCS.dataName, mode='a', complevel=4, complib='zlib')

                    print(df_expr.shape, df_EC.shape[1], df_other.shape[1], df_ttest.shape, end='\t')

            except Exception as exception:
                pass

            print(flush=True)

    if True:
        humanBatches = np.unique(ECboth[ECboth == 'Homo sapiens'].index.get_level_values('batch').values)
        mouseBatches = np.unique(ECboth[ECboth == 'Mus musculus'].index.get_level_values('batch').values)
        print(len(humanBatches), len(mouseBatches))
        
        for species in ['Homo sapiens', 'Mus musculus']:
            batches = humanBatches if species == 'Homo sapiens' else mouseBatches
            print(species, len(batches))

            fileName = 'results/DCS output/ttest_PanglaoEC.h5'

            genes = []
            genes_batches = []
            for i, batch in enumerate(os.listdir('results/DCS output/PanglaoDB/')[:]):
                if batch in batches:
                    try:
                        df_ttest = pd.read_hdf(fileName, key=batch)
                        genes.append(df_ttest.loc[df_ttest['pvalue'] <= 10**-3]['statistic'].index.values)
                        genes_batches.append(batch)
                    except Exception as exception:
                        pass

            ugenes = []
            for i, v in enumerate(genes):
                ugenes.extend(v)
            ugenes = np.unique(ugenes)
            print(len(genes), len(ugenes))

            df = pd.DataFrame(index=range(len(ugenes)), columns=genes_batches)
            for i, v in enumerate(genes):
                df.iloc[:len(v), i] = v

            df_ranks = pd.DataFrame(index=ugenes, columns=genes_batches)
            for i, v in enumerate(genes):
                df_ranks.iloc[:, i] = pd.Index(df.iloc[:, i]).reindex(df_ranks.index)[1]

            print(df_ranks)
            df_ranks.to_hdf('results/PanglaoDB_ttest_ranks_per_batch %s.h5' % species, key='df', mode='a', complevel=4, complib='zlib')

            df_ranks = df_ranks.apply(np.median, axis=1).sort_values()
            df_ranks = df_ranks[df_ranks > -1][:]
            print(df_ranks)

            df_ranks.to_excel('results/from_ranks_Panglao %s.xlsx' % species)
            