from decneo.commonFunctions import *
from decneo.analysisPipeline import Analysis
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss

def getAUCROCandPvalue(figureDataFileName, theGene):

    '''Logistic Regression and markers localization
    getAUCROCandPvalue('d:/Projects/A_Endothelial/VS/Endothelial/results/PanglaoDB_byDCS_mouse/bootstrap/All/dendrogram-heatmap-correlation-data.xlsx', 'FLT1')
    '''

    df = pd.read_excel(figureDataFileName, index_col=0, header=0)[['Markers']]
    df['Distance'] = np.abs(np.arange(len(df.index.values)) - np.where(df.index.values==theGene)[0])
    df = df.sort_values(by='Distance')

    data = df['Distance'].values[:, None]
    labels = df['Markers'].values

    AUC_ROC = roc_auc_score(~labels, data)

    clf = LogisticRegression(random_state=0).fit(data, labels)
    clf.fit(data, labels)
    alt_log_likelihood = -log_loss(labels, clf.predict_proba(data), normalize=False)

    null_prob = sum(labels) / float(labels.shape[0]) * np.ones(labels.shape)
    null_log_likelihood = -log_loss(labels, null_prob, normalize=False)
 
    p_value = scipy.stats.chi2.sf(2 * (alt_log_likelihood - null_log_likelihood), data.shape[1])
    print('AUC_ROC:', AUC_ROC, '\t', 'p-value:', p_value)

    return AUC_ROC, p_value

def hypergeometricDistributionTest(dir, file):

    '''Hypergeometric distribution function for markers in main peak using:
    scipy.stats.hypergeom(898, 22, 38).sf(10-1)
    '''

    df = pd.read_excel(dir + file, index_col=0, header=0, sheet_name='Cluster index')

    totalReceptors = df.index.values
    receptorsInPeak = getGenesOfPeak(df['Avg combo3avgs'])
    totalMarkers = df['Markers'][df['Markers'] > 0].index.values
    markersInPeak = set(receptorsInPeak).intersection(totalMarkers)

    pvalue = scipy.stats.hypergeom(len(totalReceptors), len(totalMarkers), len(receptorsInPeak)).sf(len(markersInPeak)-1)

    print('p-value:', pvalue, len(totalReceptors), len(totalMarkers), len(receptorsInPeak), len(markersInPeak))
        
    return pvalue

if __name__ == '__main__':

    '''
    Metrics comparison 1:
        All-case AUC ROC from KDR
        All-case AUC ROC from FLT1
        All-case logistic reqression p-value from KDR
        All-case logistic reqression p-value from FLT1
        All-case largest peak hypergeometric distribution test p-value
        Number of receptors with bootstrap frequency >= 0.5
        Number of receptors with bootstrap frequency >= 0.9

    Metrics comparison 2:
        Part 1: unload the bootstrap frequencies

        Part 2: compare
        Spearman correlation of bootstrap frequencies
        Pearson correlation of bootstrap frequencies
    '''

    if platform.system()=="Windows":
        wdir = 'd:/Projects/A_Endothelial/VS/Endothelial/results/' 
    else:
        wdir = '/mnt/research/piermarolab/Sergii/results/'

    if False:
        df_metricsComparison1 = pd.DataFrame()
        for majorMetric in ['spearman', 'correlation']:
            for dendrogramMetric in ['euclidean', 'correlation']:
                for linkageMethod in ['ward', 'complete', 'average']:
                    for species in ['mouse', 'human']:
                        print(1, majorMetric, dendrogramMetric, linkageMethod)

                        tempDir = wdir + 'PanglaoDB_byDCS_%s_%s/' % (species, majorMetric) + 'results majorMetric=%s, dendrogramMetric=%s, linkageMethod=%s/' % (majorMetric, dendrogramMetric, linkageMethod)
                        tempFile = 'results All %s %s %s.xlsx' % (majorMetric, dendrogramMetric, linkageMethod)
                    
                        metricsComparison1 = dict()

                        temp = getAUCROCandPvalue(tempDir + 'results All %s %s %s.xlsx' % (majorMetric, dendrogramMetric, linkageMethod), 'KDR')
                        metricsComparison1['AUCROC KDR'] = temp[0]
                        metricsComparison1['p-value KDR'] = temp[1]

                        temp = getAUCROCandPvalue(tempDir + 'results All %s %s %s.xlsx' % (majorMetric, dendrogramMetric, linkageMethod), 'FLT1')
                        metricsComparison1['AUCROC FLT1'] = temp[0]
                        metricsComparison1['p-value FLT1'] = temp[1] 

                        metricsComparison1['p-value hypergeometric'] = hypergeometricDistributionTest(tempDir, tempFile)

                        bootstrapValues = pd.read_excel(tempDir + 'Avg combo3avgs bootstrap_in-peak_genes_SD.xlsx', index_col=0, header=0)['Bootstrap'].values

                        metricsComparison1['Number of receptors >= 0.5'] = (bootstrapValues >= 0.5).sum()
                        metricsComparison1['Number of receptors >= 0.9'] = (bootstrapValues >= 0.9).sum()

                        df_metricsComparison1[(majorMetric, dendrogramMetric, linkageMethod, species)] = pd.Series(metricsComparison1)

        df_metricsComparison1.columns = pd.MultiIndex.from_tuples(df_metricsComparison1.columns.values, names=['major metric', 'dendrogram metric', 'linkage method', 'species'])
        df_metricsComparison1 = df_metricsComparison1.T
        print('\n', df_metricsComparison1)
        df_metricsComparison1.to_excel(wdir + 'df_metricsComparison1.xlsx', merge_cells=False)

    if False:
        df_metricsComparison2 = pd.DataFrame()
        seList = []
        for majorMetric in ['spearman', 'correlation']:
            for dendrogramMetric in ['euclidean', 'correlation']:
                for linkageMethod in ['ward', 'complete', 'average']:
                    for species in ['mouse', 'human']:
                        print(2, majorMetric, dendrogramMetric, linkageMethod)

                        tempDir = wdir + 'PanglaoDB_byDCS_%s_%s/' % (species, majorMetric) + 'results majorMetric=%s, dendrogramMetric=%s, linkageMethod=%s/' % (majorMetric, dendrogramMetric, linkageMethod)
                        tempFile = 'results All %s %s %s.xlsx' % (majorMetric, dendrogramMetric, linkageMethod)
                    
                        allReceptors = pd.read_excel(tempDir + tempFile, index_col=0, header=0, sheet_name='Cluster index').index

                        df_metricsComparison2[(majorMetric, dendrogramMetric, linkageMethod, species)] = pd.Series()
                        
                        seList.append(pd.read_excel(tempDir + 'Avg combo3avgs bootstrap_in-peak_genes_SD.xlsx', index_col=0, header=0)['Bootstrap'].reindex(allReceptors).fillna(0.))

        tempc = df_metricsComparison2.columns.values.copy()
        df_metricsComparison2 = pd.concat(seList, axis=1, sort=False).fillna(0.)
        df_metricsComparison2.columns = pd.MultiIndex.from_tuples(tempc, names=['major metric', 'dendrogram metric', 'linkage method', 'species'])

        df_metricsComparison2 = df_metricsComparison2.T
        print('\n', df_metricsComparison2)
        df_metricsComparison2.to_excel(wdir + 'df_metricsComparison2.xlsx', merge_cells=False)

    if False:
        df1 = pd.read_excel(wdir + 'df_metricsComparison1.xlsx', index_col=[0,1,2,3], header=0)
        df1.columns.names = ['measure']
        df1 = df1.unstack('species').reorder_levels(['species', 'measure'], axis=1).sort_index(axis=1, ascending=False)
        df1.to_excel(wdir + 'df_metricsComparison1a.xlsx', merge_cells=False)
        print(df1)

    if False:
        df2 = pd.read_excel(wdir + 'df_metricsComparison2.xlsx', index_col=[0,1,2,3], header=0).T
        df2 = df2.reorder_levels(['species', 'major metric', 'dendrogram metric', 'linkage method'], axis=1).sort_index(axis=1, ascending=False)
        
        writer = pd.ExcelWriter(wdir + 'df_metricsComparison2a.xlsx')
        df2.corr(method='pearson').to_excel(writer, 'pearson')
        df2.corr(method='spearman').to_excel(writer, 'spearman')
        writer.save()

    if True:
        df2 = pd.read_excel(wdir + 'df_metricsComparison2.xlsx', index_col=[0,1,2,3], header=0).T
        df2 = df2.reorder_levels(['species', 'major metric', 'dendrogram metric', 'linkage method'], axis=1).sort_index(axis=1, ascending=False)
        
        writer = pd.ExcelWriter(wdir + 'Metrics and 0.5 0.9 receptors lists.xlsx')
        for species in ['mouse', 'human']:
            def getFromCutoff(cutoff):

                df_temp = df2.xs(species, level='species', axis=1)

                df_temp[df_temp>=cutoff] = 1.
                df_temp[df_temp<=cutoff] = 0.

                df_temp = df_temp.loc[df_temp.sum(axis=1) > 0]

                df_temp = df_temp.sort_index()

                df_temp[(0,0,0)] = df_temp.sum(axis=1)
                df_temp = df_temp.sort_values(by=(0,0,0), ascending=False, axis=0).drop((0,0,0), axis=1)

                return df_temp

            for c in [0.9, 0.5]:
                getFromCutoff(c).to_excel(writer, '%s, cutoff %s' % (species, c))

        writer.save()