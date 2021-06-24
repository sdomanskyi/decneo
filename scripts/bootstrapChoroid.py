from decneo.commonFunctions import *
from decneo.analysisPipeline import Analysis, process

def prepareInputData_human_Choroid_remapped(sample = None):

    if platform.system() == "Windows":
        fname_all_labels = 'd:/Projects/A_Endothelial/Eye/Voigt_choroid_4567_remapped_bySD_annotation_labels_13700cells.h5'
        fname_all_expr = 'd:/Projects/A_Endothelial/Eye/Voigt_choroid_4567_remapped_bySD_DCS_all.h5'
    else:
        fname_all_labels = '/mnt/research/piermarolab/data/Eye/Voigt_choroid_4567_remapped_bySD_annotation_labels_13700cells.h5'
        fname_all_expr = '/mnt/research/piermarolab/data/Eye/Voigt_choroid_4567_remapped_bySD_DCS_all.h5'

    se_all_labels = pd.read_hdf(fname_all_labels, key='df')

    if True:
        one_index = se_all_labels[(se_all_labels == 'Endothelial')].index
        other_index = se_all_labels.index.difference(one_index)
    else:
        one_index = pd.concat([se_all_labels[(se_all_labels == 'Endothelial')].sample(n=1895, replace=False), 
                               se_all_labels[(se_all_labels == 'Macrophage')]], axis=0, sort=False).index
        other_index = se_all_labels[(se_all_labels != 'Endothelial') & 
                                    (se_all_labels != 'Macrophage') &
                                    (se_all_labels != 'Unassigned')].index

    df_all = pd.read_hdf(fname_all_expr, key='df')
    print(df_all.shape)

    df_one = df_all[one_index]

    m = (df_one.loc['LIFR']>0).astype(float).groupby(level=0).mean()
    print(m, m.mean())
 
    m = (df_one.loc['KDR']>0).astype(float).groupby(level=0).mean()
    print(m, m.mean())   

    m = (df_one.loc['PECAM1']>0).astype(float).groupby(level=0).mean()
    print(m, m.mean())

    exit()

    if not sample is None:
        df_one = df_one.sample(n=sample, axis=1, replace=False)

    df_one = df_one.loc[df_one.sum(axis=1) > 0.]

    print(df_one)
    print(df_one.columns.shape[0], dict(zip(*np.unique(df_one.columns.get_level_values('batch'), return_counts=True))))
    print(df_one.columns.shape[0], dict(zip(*np.unique(se_all_labels[df_one.columns].values, return_counts=True))))

    df_other = df_all[other_index].reindex(df_one.index)
    
    if not sample is None:
        df_other = df_other.sample(n=sample, axis=1, replace=False)

    print(df_other)
    print(df_other.columns.shape[0], dict(zip(*np.unique(se_all_labels[df_other.columns].values, return_counts=True))))
    print(df_other.columns.shape[0], dict(zip(*np.unique(df_other.columns.get_level_values('batch'), return_counts=True))))

    #df_one.to_hdf('demo/demoChoroidData.h5', key='dfa', mode='a', complevel=4, complib='zlib')
    #df_other.to_hdf('demo/demoChoroidData.h5', key='dfb', mode='a', complevel=4, complib='zlib')

    return df_one, df_other

if __name__ == '__main__':

    if platform.system()=="Windows":
        wdir = 'd:/Projects/A_Endothelial/VS/Endothelial/results/' 
    else:
        wdir = '/mnt/research/piermarolab/Sergii/results/'

    #prepareInputData_human_Choroid_remapped()

    #process(*prepareInputData_human_Choroid_remapped(), *(None, None),
    #        wdir + 'choroid Voigt remapped test 10/', wdir + 'PanglaoDB_byDCS_mouse/bootstrap/All/', 
    #        nCPUs=4 if platform.system()=="Windows" else 8, parallelBootstrap=True,
    #        genesOfInterest=receptorsListHugo_2555, knownRegulators=gEC22, exprCutoff1=0.01, perEachOtherCase=False,
    #        nBootstrap=10, dendrogramMetric = 'euclidean', dendrogramLinkageMethod = 'average')

    if False:
        process(*(None, None), *(None, None),
                wdir + 'choroid Voigt remapped/', wdir + 'PanglaoDB_byDCS_mouse/bootstrap/All/', 
                nCPUs=4 if platform.system()=="Windows" else 8, parallelBootstrap=True,
                genesOfInterest=receptorsListHugo_2555, knownRegulators=gEC22, exprCutoff1=0.01, perEachOtherCase=False,
                nBootstrap=10, majorMetric='spearman', dendrogramMetric = 'euclidean', dendrogramLinkageMethod = 'average')
        
    if True:


        Analysis(workingDir=wdir + 'choroid Voigt remapped/', otherCaseDir=wdir + 'PanglaoDB_byDCS_mouse/bootstrap/All/', genesOfInterest=receptorsListHugo_2555, knownRegulators=gEC22, panels = ['fraction', 'binomial', 'top50', 'combo3avgs', 'max'], externalPanelsData=externalPanelsData).reanalyzeMain(togglePublicationFigure=True, includeClusterNumber=False, toggleIncludeHeatmap=False, toggleCalculateMeasures=False, toggleExportFigureData=False)

    if False:
        np.random.seed(0)
        m, n, k = 15000, 5000, 500
        df_sample = pd.DataFrame(index=['Gene %s' % i for i in range(m)], columns=['Cell %s' % i for i in range(n)], data=np.random.rand(m, n))
        df_sample.iloc[0] += 3. * df_sample.iloc[1]
        df_sample.iloc[2] += 5. * df_sample.iloc[-10]
        print(df_sample.shape)

        temp_genes = pd.Index(['Gene %s' % i for i in range(k)])

        #input('*')
        sT = time.time()
        measure = cdist(df_sample.loc[temp_genes].apply(lambda s: pd.Series(index=s.index, data=scipy.stats.rankdata(s.values)), axis=1),
                       df_sample.apply(lambda s: pd.Series(index=s.index, data=scipy.stats.rankdata(s.values)), axis=1),
                       metric='correlation').T
        df_measure = pd.DataFrame(data=measure, index=df_sample.index, columns=temp_genes)
        print(df_measure.round(6))
        print('**', time.time() - sT)

        #input('*')
        sT = time.time()
        measure = 1. - scipy.stats.spearmanr(df_sample.loc[temp_genes].values, df_sample.values, axis=1).correlation[-len(df_sample.values):, :len(temp_genes)]
        df_measureL = pd.DataFrame(data=measure, index=df_sample.index, columns=temp_genes)
        print(df_measureL.round(6))
        print('**', time.time() - sT)

        print('isGood:', ((df_measure.values - df_measureL.values) < 10**-8).all())

        #input('**')
        sT = time.time()
        measure = cdist(df_sample.loc[temp_genes].values, df_sample.values, metric='correlation').T
        df_measure = pd.DataFrame(data=measure, index=df_sample.index, columns=temp_genes)
        print(df_measure)
        print('***', time.time() - sT)
