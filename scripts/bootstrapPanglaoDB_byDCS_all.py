from commonFunctions import *
from analysisPipeline import Analysis

def preparePerSampleDEgenes():

    wdir = '/mnt/research/piermarolab/Sergii/DCS output/'

    processedCollectedEC = wdir + 'PanglaoDB_EC_byDCS_normed_and_filtered.h5'
    tempECindex = wdir + 'PanglaoDB_ECcells.h5'
    fileNameTtest = wdir + 'ttest_PanglaoEC.h5'
    fileNameEC = wdir + 'PanglaoDB_EC.h5'
    fileNameOther = wdir + 'PanglaoDB_nonEC.h5'
    PanglaoDB_DCS_dir = wdir + 'PanglaoDB/'

    # Prepare index of EC
    if True:
        ECcells1 = pd.Series(index=pd.read_hdf(processedCollectedEC, key='Homo sapiens').columns, data='Homo sapiens')
        ECcells2 = pd.Series(index=pd.read_hdf(processedCollectedEC, key='Mus musculus').columns, data='Mus musculus')
        ECboth = pd.concat([ECcells1, ECcells2], axis=0, sort=False)
        ECboth.to_hdf(tempECindex, key='df', mode='a', complevel=4, complib='zlib')
        print(ECboth)

        ECboth = pd.read_hdf(tempECindex, key='df')

    # Calculate ttest for EC and non-EC populations for each SRS (1368 of them). Also save expression of the cells used
    if True:
        import DigitalCellSorter

        for i, batch in enumerate(os.listdir(PanglaoDB_DCS_dir)):
            print('Processing', i, batch, end='\t')
            DCS = DigitalCellSorter.DigitalCellSorter(dataName=batch, saveDir=PanglaoDB_DCS_dir + batch)

            try:
                ECcells = ECboth.xs(key=DCS.dataName, axis=0, level='batch', drop_level=False).index

                if len(ECcells) >= 10:
                    DCS.loadExpressionData()
                    df_expr = DCS.df_expr.xs(key=DCS.dataName, axis=1, level='batch', drop_level=False)
                    columns = pd.MultiIndex.from_arrays([df_expr.columns.get_level_values('batch'), df_expr.columns.get_level_values('cell')])
                    df_expr = pd.DataFrame(data=df_expr.values, index=df_expr.index, columns=columns)

                    df_EC = df_expr[ECcells]
                    df_other = df_expr[df_expr.columns.difference(df_EC.columns)]
                    if df_other.shape[1] > 1000:
                        df_other = df_other.sample(n=1000, axis=1)

                    df_EC.to_hdf(fileNameEC, key=DCS.dataName, mode='a', complevel=4, complib='zlib')
                    df_other.to_hdf(fileNameOther, key=DCS.dataName, mode='a', complevel=4, complib='zlib')
                
                    df_ttest = pd.DataFrame(index=df_EC.index, columns=['statistic', 'pvalue'])
                    ttest = scipy.stats.ttest_ind(df_EC.values, df_other.values, axis=1)
                    df_ttest['statistic'] = ttest[0]
                    df_ttest['pvalue'] = ttest[1]
                    df_ttest = df_ttest.sort_values('statistic', ascending=False).dropna()
                    df_ttest.to_hdf(fileNameTtest, key=DCS.dataName, mode='a', complevel=4, complib='zlib')

                    print(df_expr.shape, df_EC.shape[1], df_other.shape[1], df_ttest.shape, end='\t')

            except Exception as exception:
                pass

            print(flush=True)

    # Determine DE gene ranks
    if True:
        humanBatches = np.unique(ECboth[ECboth == 'Homo sapiens'].index.get_level_values('batch').values)
        mouseBatches = np.unique(ECboth[ECboth == 'Mus musculus'].index.get_level_values('batch').values)
        print(len(humanBatches), len(mouseBatches))
        
        for species in ['Homo sapiens', 'Mus musculus']:
            batches = humanBatches if species == 'Homo sapiens' else mouseBatches
            print(species, len(batches))

            genes = []
            genes_batches = []
            for i, batch in enumerate(os.listdir(PanglaoDB_DCS_dir)):
                if batch in batches:
                    try:
                        df_ttest = pd.read_hdf(fileNameTtest, key=batch)
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
            df_ranks.to_hdf('/mnt/research/piermarolab/Sergii/DCS output/PanglaoDB_DCS_ttest_ranks_per_batch %s.h5' % species, key='df', mode='a', complevel=4, complib='zlib')

            df_ranks = df_ranks.apply(np.median, axis=1).sort_values()
            df_ranks = df_ranks[df_ranks > -1][:]
            print(df_ranks)

    return

def prepareInput(fileName, species):

    df = pd.read_hdf('/mnt/home/domansk6/Projects/Endothelial/results/DCS output/PanglaoDB_EC_byDCS_normed_and_filtered.h5', key=species)
    df.to_hdf(fileName, key='df', mode='a', complevel=4, complib='zlib')
    print(df)
    
    df_ranks = pd.read_hdf('/mnt/research/piermarolab/Sergii/DCS output/PanglaoDB_DCS_ttest_ranks_per_batch %s.h5' % species, key='df')
    df_ranks.to_hdf(fileName, key='df_ranks', mode='a', complevel=4, complib='zlib')
    print(df_ranks)

    return

if __name__ == '__main__':

    if False:
        preparePerSampleDEgenes()

    args = dict(genesOfInterest=receptorsListHugo_2555, knownRegulators=gEC23, nCPUs=4 if platform.system()=="Windows" else 20, 
                panels=['combo3avgs', 'combo4avgs', 'fraction', 'binomial', 'markers', 'top50'], nBootstrap=100, perEachOtherCase=True)

    dirHuman = '/mnt/home/domansk6/Projects/Endothelial/results/PanglaoDB_byDCS_human/'
    dirMouse = '/mnt/home/domansk6/Projects/Endothelial/results/PanglaoDB_byDCS_mouse/'

    anHuman = Analysis(**dict(args, workingDir=dirHuman, otherCaseDir=dirMouse))
    prepareInput(anHuman.dataSaveName, 'Homo sapiens')
    anHuman.preparePerBatchCase(exprCutoff=0.05)
    anHuman.prepareBootstrapExperiments()
    anHuman.analyzeBootstrapExperiments()

    anMouse = Analysis(**dict(args, workingDir=dirMouse, otherCaseDir=dirHuman))
    prepareInput(anMouse.dataSaveName, 'Mus musculus')
    anMouse.preparePerBatchCase(exprCutoff=0.05)
    anMouse.prepareBootstrapExperiments()
    anMouse.analyzeBootstrapExperiments()

    anHuman.analyzeBootstrapExperiments()

    anHuman.reanalyzeMain()
    anHuman.analyzeCombinationVariant('Avg combo3avgs')
    anHuman.analyzeCombinationVariant('Avg combo4avgs')
    anHuman.analyzeAllPeaksOfCombinationVariant('Avg combo3avgs', nG=8, nE=30, fcutoff=0.5, width=50)
    anHuman.analyzeAllPeaksOfCombinationVariant('Avg combo4avgs', nG=8, nE=30, fcutoff=0.5, width=50)
    anHuman.scramble(['Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo3/', M=20)  
    anHuman.scramble(['Markers', 'Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo4/', M=20)  

    anMouse.reanalyzeMain()
    anMouse.analyzeCombinationVariant('Avg combo3avgs')
    anMouse.analyzeCombinationVariant('Avg combo4avgs')
    anMouse.analyzeAllPeaksOfCombinationVariant('Avg combo3avgs', nG=8, nE=30, fcutoff=0.5, width=50)
    anMouse.analyzeAllPeaksOfCombinationVariant('Avg combo4avgs', nG=8, nE=30, fcutoff=0.5, width=50)
    anMouse.scramble(['Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo3/', M=20)  
    anMouse.scramble(['Markers', 'Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo4/', M=20)  