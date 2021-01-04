from decneo.commonFunctions import *
from decneo.analysisPipeline import Analysis

if __name__ == '__main__':

    # Checking number og GO terms and terms-proteins annotations
    if False:
        import json 
        with gzip.GzipFile('d:/Projects/A_Endothelial/VS/Endothelial/dev/geneLists/humanGeneOntAssoc.json.gz', 'r') as tempFile:
            inputFile = json.loads(tempFile.read().decode('utf-8'))[1]

        with gzip.open('d:/Projects/A_Endothelial/VS/Endothelial/dev/geneLists/goBasicObo.txt.gz', 'r') as tempFile:
            inputFile = tempFile.read().decode().split('id: GO:')
            terms = dict()
            terms.update(dict(biological_process = [item for item in inputFile if 'namespace: biological_process' in item]))
            terms.update(dict(cellular_component = [item for item in inputFile if 'namespace: cellular_component' in item]))
            terms.update(dict(molecular_function = [item for item in inputFile if 'namespace: molecular_function' in item]))

            total = 0
            for key in terms.keys():
                temp = len(terms[key])
                print(key, '\t', temp, 'terms')
                total += temp

            print('Total terms:\t', total)

    # II. Re-calculate dendrogram randomization statistics
    if False:
        measures3 = ['Binomial -log(pvalue)', 'Top50 overlap', 'Fraction']

        if platform.system() == "Windows":
            nCPUs = 2 
            dir = 'results/PanglaoDB_byDCS_mouse/'
        else:
            nCPUs = 100
            dir = '/mnt/research/piermarolab/Sergii/results/PanglaoDB_byDCS_mouse/'

        try:
            ind = sys.argv[1]
        except Exception as exception:
            ind = 0

        import time
        sT = time.time()
        Analysis(dir, dir + '', nCPUs=nCPUs).scramble(measures3, case='All', subDir='All/combo3/', N=10**2, M=10**6, getMax=True, maxSuff='x1000000_%s' % ind)
        print(np.round(time.time() - sT, 3))
 
    # Re-calculate dendrogram randomization statistics
    if False:
        def run(dir, case):
            print(dir, case)

            an = Analysis(dir, dir + '', nCPUs=20)

            an.scramble(['Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], case=case, subDir='%s/combo3/' % case, M=100)

            return

        try:
            case = sys.argv[1]
        except Exception as exception:
            case = 'All'

        run('/mnt/research/piermarolab/Sergii/results/PanglaoDB_byDCS_mouse/', case)

        #run('/mnt/research/piermarolab/Sergii/results/PanglaoDB_byDCS_mouse/')
        #run('/mnt/research/piermarolab/Sergii/results/PanglaoDB_byDCS_human/')
        #run('/mnt/research/piermarolab/Sergii/PanglaoDB_byAlona/PanglaoDB_byDCS_mouse/')
        #run('/mnt/research/piermarolab/Sergii/PanglaoDB_byAlona/PanglaoDB_byDCS_human/')

        pass

    # Overlay groups on main case
    if False:
        print(os.getcwd())
        se_measure = pd.read_excel('dev/DCS mouse all combo3.xlsx', sheet_name='Sheet1', index_col=0, header=0)['Avg combo3avgs']
        print(se_measure)

        se_groups = pd.read_excel('dev/DCS mouse all combo3.xlsx', sheet_name='Sheet2', index_col=0, header=0)['Genes'].str.split(', ').apply(pd.Series)
        se_groups.index = se_groups.index.str.replace('P', '').astype(int)
        se_groups = se_groups.stack().reset_index().set_index(0)['Group'].reindex(se_measure.index).fillna(0).astype(int)
        print(se_groups)

        fig, ax = plt.subplots()

        ax.plot(range(len(se_measure.index)), se_measure.values)

        u = np.unique(se_groups.values)
        for group in u:
            if group > 0:
                x = np.where(se_groups.values==group)[0]

                c = plt.cm.jet(group/len(u))

                label = 'P%s' % group

                ax.bar(x, se_measure.values[x], width=1.0, color=c, label=label)

                location = int(np.median(x)) 
                ax.text(location, se_measure.values[location], label, fontsize=14, va='center', ha='center').set_path_effects([path_effects.Stroke(linewidth=2., foreground='w'), path_effects.Normal()])

        ax.set_xticks([])
        ax.set_xticklabels([])
        yticks = (0, 0.03)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks)
        ax.tick_params(axis='y', labelsize=16, width=0.75, length=3)
        ax.text(-0.05, 0.5, 'Combination of 3', fontsize=16, rotation=90, va='center', ha='right', transform=ax.transAxes)
        ax.set_xlim([0, len(se_measure)])
        ax.set_ylim([0, 0.03])

        plt.legend()

        plt.show()
        plt.close(fig)

    # Collect Supplementary Table 6
    if False:
        wdir = 'd:/Projects/A_Endothelial/VS/Endothelial/results/PanglaoDB_byDCS_mouse/'

        writer = pd.ExcelWriter(wdir + 'Supplementary Table 6.xlsx')

        for file in ['GO11.xlsx', 'KEGG11.xlsx', 'Reactome11.xlsx']:
            print(file)

            dfs = []
            for sheet in range(1,11+1):
                print(sheet)
                df_temp = pd.concat([pd.read_excel(wdir + file, sheet_name=str(sheet), index_col=0, header=0)], keys=['P%s' % sheet], names=['Group'], sort=False)
                dfs.append(df_temp)

            print('Combining')
            dfs = pd.concat(dfs, axis=0, sort=False)
            print(dfs)

            dfs.to_excel(writer, file.split('11.xlsx')[0], merge_cells=False)

        writer.save()

    # Markers localization p-value calculation
    if False:
        df_in = pd.read_excel('d:/Projects/A_Endothelial/VS/Endothelial/results/PanglaoDB_byDCS_mouse_correlation/bootstrap/All/dendrogram-heatmap-correlation-data.xlsx', index_col=0, header=0)[['Markers']]

        print('Number of receptors:', len(df_in))

        for gene in ['KDR', 'FLT1']:
            print('\n\n', gene, ':')

            df = df_in.copy()
            df['Distance'] = np.abs(np.arange(len(df.index.values)) - np.where(df.index.values==gene)[0])

            data = df['Distance'].values[:, None]
            labels = df['Markers'].values

            np.random.seed(0)
            np.random.shuffle(labels)

            if True:
                from rpy2.robjects import r as R
                from rpy2.robjects import IntVector, FloatVector, DataFrame

                try:
                    from rpy2.robjects.packages import importr
                    utils = importr('dplyr')
                    utils = importr('gravity')
                except Exception as exception:
                    print(exception)
                    import rpy2.robjects.packages as rpackages
                    utils = rpackages.importr('utils')
                    utils = rpackages.importr('gravity')
                    utils.chooseCRANmirror(ind=1)
                    utils.install_packages('dplyr')
                    utils.install_packages('gravity')
                finally:
                    from rpy2.robjects.packages import importr
                    utils = importr('dplyr')
                    utils = importr('gravity')

                dataf = DataFrame({'label': IntVector(tuple(labels)), 'distance': FloatVector(tuple(data.T[0]))})
                fit = R('function(x) ppml(dependent_variable="label", distance="distance", additional_regressors=NULL, robust=TRUE, data=x)')(dataf)
                #print(fit)

                # Deviance is -2.*log_likelihood
                altDeviance = list(fit[9].items())[0][1]
                nullDeviance = list(fit[11].items())[0][1]
                p_value = scipy.stats.chi2.sf(nullDeviance - altDeviance, 1)
                print('Poisson pseudo-maximum likelihood estimation (PPML) by J.M.C. Santos Silva & Silvana Tenreyro, 2006.')
                print('Implemented in R in "gravity: Estimation Methods for Gravity Models" at:')
                print('https://rdrr.io/cran/gravity/man/ppml.html')
                print('Robust PPML based (a.k.a. QMLE) deviances and their difference test (chi2 p-value):\n\t', 
                      'Null deviance:\t', np.round(nullDeviance, 1), '\n\t',
                      'Alternative deviance:\t', np.round(altDeviance, 1), '\n\t',
                      'p-value:\t', '%.1e' % p_value, '\n')

            if True:
                from sklearn.linear_model import LogisticRegression
                from sklearn.metrics import log_loss
                clf = LogisticRegression(random_state=0).fit(data, labels)
                clf.fit(data, labels)
                prprob = clf.predict_proba(data)
                alt_log_likelihood = -log_loss(labels, prprob, normalize=False)
                null_prob = sum(labels) / float(labels.shape[0]) * np.ones(labels.shape)
                null_log_likelihood = -log_loss(labels, null_prob, normalize=False)
                altDevianceLogit = -2.*alt_log_likelihood
                nullDevianceLogit = -2.*null_log_likelihood
                p_value_logit = scipy.stats.chi2.sf(nullDevianceLogit - altDevianceLogit, 1)
                print('Previously I used:')
                print('Logistic regression based Deviances and their difference test (chi2 p-value):\n\t', 
                      'Null deviance:\t', np.round(nullDevianceLogit, 1), '\n\t',
                      'Alternative deviance:\t', np.round(altDevianceLogit, 1), '\n\t',
                      'p-value:\t', '%.1e' % p_value_logit, '\n')

                AUC_ROC = roc_auc_score(~labels, data)
                print('AUC_ROC:', np.round(AUC_ROC, 2))

            if False:
                from statsmodels.graphics.gofplots import qqplot
                qqplot(prprob.T[0])
                plt.show()

    # Calculate Table of correlation comparing methods
    if False:
        #wdir = 'd:/Projects/A_Endothelial/VS/Endothelial/results/'
        #anMouse = Analysis(wdir + 'PanglaoDB_byDCS_mouse/', wdir + '', genesOfInterest=receptorsListHugo_2555, knownRegulators=gEC22)
        #anMouse.analyzeAllPeaksOfCombinationVariant('Avg combo3avgs', nG=15, nE=15, fcutoff=0.5, width=50)

        df = pd.read_excel('d:/Projects/A_Endothelial/VS/Endothelial/results/for meeting 10 15 2020/Table 1 for paper - for corr.xlsx', sheet_name='Table 1', index_col=0, header=[0,1,2])
        print(df.shape)

        df[df < 0.5] = 0.
        df_corr = df.corr()
        np.fill_diagonal(df_corr.values, np.nan)
        df_corr.to_excel('d:/Projects/A_Endothelial/VS/Endothelial/results/for meeting 10 15 2020/Corr of Table 1.xlsx')
        print(df_corr)

    # Hypergeometric test
    if True:
        # Markers in main peak
        print(scipy.stats.hypergeom(898, 22, 38).sf(10-1)) # 2.4e-09
        
        # Markers in main bootstrap frequency peak
        print(scipy.stats.hypergeom(898, 22, 18).sf(7-1))  # 5.0e-08   
        
        # Overlap in split by SRA
        print(scipy.stats.hypergeom(898, 21, 23).sf(17-1))  # 1.5e-27   
        
        # Overlap in split by SRS
        print(scipy.stats.hypergeom(898, 26, 20).sf(16-1))  # 3.3e-24   
        
        # Overlap in 2nd validation and main
        print(scipy.stats.hypergeom(1099, 18, 25).sf(16-1))  # 1.6e-27