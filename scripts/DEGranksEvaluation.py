from decneo.commonFunctions import *
from decneo.analysisPipeline import process

def getDEG(dfa, dfb, method='ttest'):

    '''
    method = {'ttest', 'wilcoxon', 'mannwhitneyu'}
    '''

    genes = []
    batches = []
    for batch in np.unique(dfa.columns.get_level_values('batch').values):

        df_temp_a = dfa.xs(batch, level='batch', axis=1, drop_level=False)
        df_temp_b = dfb.xs(batch, level='batch', axis=1, drop_level=False)

        if method == 'ttest':
            ttest = scipy.stats.ttest_ind(df_temp_a.values, df_temp_b.values, axis=1)
            df_test = pd.DataFrame(index=dfa.index, columns=['statistic', 'pvalue'])
            df_test['statistic'] = ttest[0]
            df_test['pvalue'] = ttest[1]

        elif method == 'wilcoxon':
            min_size = min(df_temp_a.shape[1], df_temp_b.shape[1])
            df_temp_a = df_temp_a.sample(min_size, axis=1)
            df_temp_b = df_temp_b.sample(min_size, axis=1)
            df_temp = pd.concat([df_temp_a, df_temp_b], axis=1, sort=False)
            df_temp = df_temp.loc[df_temp.sum(axis=1) > 0]
            df_test = df_temp.apply(lambda v: pd.Series(np.array(scipy.stats.wilcoxon(v[:min_size], v[min_size:]))), axis=1)
            df_test.columns = ['statistic', 'pvalue']

        elif method == 'mannwhitneyu':
            df_temp_a = df_temp_a.apply(np.array, axis=1)
            df_temp_b = df_temp_b.apply(np.array, axis=1)
            df_temp = pd.concat([df_temp_a, df_temp_b], axis=1, sort=False)   
            df_temp = df_temp.loc[(df_temp[0].apply(np.sum) + df_temp[1].apply(np.sum)) > 0]
            df_test = df_temp.apply(lambda v: pd.Series(np.array(scipy.stats.mannwhitneyu(v[0], v[1]))), axis=1)
            df_test.columns = ['statistic', 'pvalue']

        #elif method == 'mannwhitneyu':
        #    df_temp_a = df_temp_a.apply(np.array, axis=1)
        #    df_temp_b = df_temp_b.apply(np.array, axis=1)
        #    df_temp = pd.concat([df_temp_a, df_temp_b], axis=1, sort=False)   
        #    df_temp = df_temp.loc[(df_temp[0].apply(np.sum) + df_temp[1].apply(np.sum)) > 0]

        #    def mfunc(u, v):

        #        u, v = u[u>0], v[v>0]

        #        if len(u) >= 10 and len(v) >= 10:
        #            return scipy.stats.mannwhitneyu(u, v)
        #        else:
        #            return 0., 1.

        #    df_test = df_temp.apply(lambda v: pd.Series(np.array(mfunc(v[0], v[1]))), axis=1)
        #    df_test.columns = ['statistic', 'pvalue']

        df_ttest = df_ttest.sort_values('statistic', ascending=False).dropna()
        genes.append(df_ttest.loc[df_ttest['pvalue'] <= 10**-3]['statistic'].index.values)
        batches.append(batch)

    ugenes = []
    for i, v in enumerate(genes):
        ugenes.extend(v)
    ugenes = np.unique(ugenes)

    df = pd.DataFrame(index=range(len(ugenes)), columns=batches)
    for i, v in enumerate(genes):
        df.iloc[:len(v), i] = v

    df_ranks = pd.DataFrame(index=ugenes, columns=batches)
    for i, v in enumerate(genes):
        df_ranks.iloc[:, i] = pd.Index(df.iloc[:, i]).reindex(df_ranks.index)[1]

    return df_ranks

if __name__ == '__main__':
    
    wdir = 'd:/Projects/A_Endothelial/VS/Endothelial/results/'
    demoData = 'd:/Projects/A_Endothelial/VS/Endothelial/data/VoightChoroid4567RemappedData.h5'
    ranksFile = 'd:/Projects/A_Endothelial/VS/Endothelial/results/df_ranks_ttest_and_wilcoxon.h5'

    # Calculate DEG ranks by 3 methods
    if False:
        dfs = pd.read_hdf(demoData, key='dfa'), pd.read_hdf(demoData, key='dfb')
    
        for method in ['ttest', 'wilcoxon', 'mannwhitneyu']:
            df_ranks = getDEG(*dfs, method=method)
            df_ranks.to_hdf(ranksFile, key=method, **phdf)
            print(df_ranks)

    # Testing demo with 3 DEG methods
    if False:
        method = 'mannwhitneyu' # trying each one of {'ttest', 'mannwhitneyu'}

        process(pd.read_hdf(demoData, key='dfa'),   # Endothelial cells
                pd.read_hdf(demoData, key='dfb'),   # Non-endothelial cells
                None, None,                         # Comparison dataset is provided
                wdir + 'demo/%s_3/' % method,               
                wdir + 'demo/fromPanglaoDBmouseAllbyDCS/',
                nBootstrap=1,
                nCPUs = 2,
                methodForDEG=method,                
                parallelBootstrap=True,             # Set False if RAM is limited
                exprCutoff1=0.01,                   # Gene expression cutoff
                perEachOtherCase=False)

    # Data normality testing
    if False:
        for key in ['dfa', 'dfb']:
            df = pd.read_hdf(demoData, key=key)
            df = df.loc[(df>0).sum(axis=1) >= 20].apply(lambda s: scipy.stats.normaltest(s.values[s.values > 0])[1], axis=1)
            print(df, '\n', df.mean(), '\n')

    if True:
        cols = ['ttest', 'mannwhitneyu_3']
        ars = [pd.read_hdf(wdir + 'demo/%s/data.h5' % subDir, key='df_ranks').apply(np.median, axis=1).sort_values(ascending=False)[:1000].sort_index().index.values for subDir in cols]

        dfal = pd.DataFrame(index=np.unique(np.hstack(ars)), columns=cols)
        for i, ar in enumerate(ars):
            dfal[dfal.columns[i]] = pd.Series(index=ar, data=1).reindex(dfal.index)

        print(dfal)
        dfal.to_excel(wdir + 'dfal2.xlsx')

