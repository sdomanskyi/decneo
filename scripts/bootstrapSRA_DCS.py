from scRegulation.commonFunctions import *
from scRegulation.analysisPipeline import Analysis

if __name__ == '__main__':

    if False:
        cases = pd.read_excel('DCS tissues SRAs.xlsx', sheet_name='Selected DCS', index_col=[0,1,2], header=0)
        cases = cases[cases['cells'] >= 100]
        cases = cases['SRS'].str.replace(' ', '').str.split(',')
        print(cases)

        for (tissue, species, SRA), SRSs in cases.iteritems():
            print('\n', tissue, species, SRA, len(SRSs), flush=True)

            args = dict(genesOfInterest=receptorsListHugo_2555, knownRegulators=gEC23, nCPUs=4 if platform.system()=="Windows" else 10, panels=['combo3avgs', 'combo4avgs', 'fraction', 'binomial', 'markers', 'top50'], nBootstrap=100, perEachOtherCase=False)

            allPanglaoDBmouse = '/mnt/home/domansk6/Projects/Endothelial/results/PanglaoDB_byDCS_mouse/bootstrap/All/'
            allPanglaoDBhuman = '/mnt/home/domansk6/Projects/Endothelial/results/PanglaoDB_byDCS_human/bootstrap/All/'
        
            an = Analysis(**dict(args, workingDir='/mnt/research/piermarolab/Sergii/DCS_SRAs/PanglaoDB_byDCS %s %s %s/' % (tissue, species, SRA), otherCaseDir=allPanglaoDBmouse if species == 'Homo sapiens' else allPanglaoDBhuman))

            if False:
                if os.path.isfile(an.workingDir + 'results.png'):
                    continue

                if not os.path.isfile(an.dataSaveName):
                    df_EC = pd.concat([pd.read_hdf('/mnt/research/piermarolab/Sergii/DCS output/PanglaoDB_EC.h5', key=SRA + '_' + SRS) for SRS in SRSs], axis=1, sort=False).fillna(0.)
                    df_other = pd.concat([pd.read_hdf('/mnt/research/piermarolab/Sergii/DCS output/PanglaoDB_nonEC.h5', key=SRA + '_' + SRS) for SRS in SRSs], axis=1, sort=False).fillna(0.)

                    if len(SRSs) < 5:
                        nBatches = 10
                        df_EC.columns = pd.MultiIndex.from_arrays([np.random.permutation(np.hstack([np.array([str(i)]*len(v), dtype=str) for i, v in enumerate(np.array_split(df_EC.columns.get_level_values('batch').values, nBatches))])), df_EC.columns.get_level_values('cell')], names=['batch', 'cell'])
                        df_other.columns = pd.MultiIndex.from_arrays([np.random.permutation(np.hstack([np.array([str(i)]*len(v), dtype=str) for i, v in enumerate(np.array_split(df_other.columns.get_level_values('batch').values, nBatches))])), df_other.columns.get_level_values('cell')], names=['batch', 'cell'])

                    print(df_EC, df_other)
                    df_EC = df_EC.T.loc[~df_EC.T.index.duplicated(keep='first')].T
                    df_other = df_other.T.loc[~df_other.T.index.duplicated(keep='first')].T
                    print(df_EC.shape, df_other.shape)

                    an.prepareDEG(df_EC, df_other)

                an.preparePerBatchCase(exprCutoff=0.05)
                an.prepareBootstrapExperiments()
                an.analyzeBootstrapExperiments()
                an.reanalyzeMain()
                an.analyzeCombinationVariant('Avg combo3avgs')
                an.analyzeAllPeaksOfCombinationVariant('Avg combo3avgs', nG=8, nE=30, fcutoff=0.5, width=50)
                an.analyzeCombinationVariant('Avg combo4avgs')
                an.analyzeAllPeaksOfCombinationVariant('Avg combo4avgs', nG=8, nE=30, fcutoff=0.5, width=50)
                an.scramble(['Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo3/', M=20)  
                an.scramble(['Markers', 'Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo4/', M=20)  

            else:
                name = '%s %s %s' % (tissue, species, SRA)

                # Collect bootstrap counts for combo3 and combo4
                b3 = pd.read_excel(an.workingDir + 'Avg combo3avgs_variant.xlsx', index_col=0)['Bootstrap']
                b4 = pd.read_excel(an.workingDir + 'Avg combo4avgs_variant.xlsx', index_col=0)['Bootstrap']
                b3.name = 'combo3'
                b4.name = 'combo4'
                b3.to_hdf('results_DCS_35_SRA.h5', key='combo3/b/%s' % name, mode='a', complevel=4, complib='zlib')
                b4.to_hdf('results_DCS_35_SRA.h5', key='combo4/b/%s' % name, mode='a', complevel=4, complib='zlib')

                # Collect max random values for combo3 and combo4
                r3 = pd.read_excel(an.workingDir + 'random/combo3/se_distribution.xlsx', index_col=0)[0]
                r4 = pd.read_excel(an.workingDir + 'random/combo4/se_distribution.xlsx', index_col=0)[0]
                r3.name = 'combo3'
                r4.name = 'combo4'
                r3.to_hdf('results_DCS_35_SRA.h5', key='combo3/r/%s' % name, mode='a', complevel=4, complib='zlib')
                r4.to_hdf('results_DCS_35_SRA.h5', key='combo4/r/%s' % name, mode='a', complevel=4, complib='zlib')

                # Collect sub-peaks lists for combo3 and combo4
                s3 = pd.read_excel(an.workingDir + 'All peaks Avg combo3avgs. nG8-nE30.xlsx', index_col=0)['genes'].str.replace(' ', '').str.split(',')
                s4 = pd.read_excel(an.workingDir + 'All peaks Avg combo4avgs. nG8-nE30.xlsx', index_col=0)['genes'].str.replace(' ', '').str.split(',')
                s3.name = 'combo3'
                s3.to_hdf('results_DCS_35_SRA.h5', key='combo3/s/%s' % name, mode='a', complevel=4, complib='zlib')
                s4.to_hdf('results_DCS_35_SRA.h5', key='combo4/s/%s' % name, mode='a', complevel=4, complib='zlib')

    if True:
        def getDf(keys, combo, parameter, species):

            dfs = []
            for name in np.unique(keys.T[2]):
                if species in name:
                    b = pd.read_hdf('results_DCS_35_SRA.h5', key='/%s/%s/%s' % (combo, parameter, name))
                    b.name = name
                    dfs.append(b)

            dfs = pd.concat(dfs, axis=1, sort=False).fillna(0.)
            print(dfs.shape)

            return dfs

        def plotDf(df, figName):
            fig = plt.figure(figsize=(8,10))

            n_clusters = 10
            ax = fig.add_axes([0.35, 0.7, 0.55, 0.15], frame_on=False)
            Z = hierarchy.linkage(np.nan_to_num(df.values, nan=df.values.max()), method='ward', optimal_ordering=True)
            origLineWidth = matplotlib.rcParams['lines.linewidth']
            matplotlib.rcParams['lines.linewidth'] = 0.5
            D = hierarchy.dendrogram(Z, ax=ax, color_threshold=0, above_threshold_color='k', orientation='top')
            matplotlib.rcParams['lines.linewidth'] = origLineWidth
            ax.set_xticklabels([])
            ax.set_xticks([])
            ax.set_yticklabels([])
            ax.set_yticks([])

            plt.title(figName)

            ax = fig.add_axes([0.35, 0.3, 0.55, 0.4], frame_on=False)
            cmap = plt.cm.hot
            cmap.set_bad('grey')
            im = ax.imshow(np.ma.array(df.values[:,D['leaves']], mask=np.isnan(df.values)), cmap=cmap, aspect='auto', interpolation='None', extent=(-0.5, df.shape[0] - 0.5, df.shape[1] - 0.5, -0.5))
            ax.set_xticks(range(len(df.columns)))
            ax.set_yticks(range(len(df.columns)))
            ax.set_xticklabels(df.columns.values[D['leaves']])
            ax.set_yticklabels(df.columns.values)
            ax.tick_params(axis='x', labelsize=8, width=0.25, length=1, rotation=90)
            ax.tick_params(axis='y', labelsize=8, width=0.25, length=1, rotation=0)

            ax = fig.add_axes([0.9, 0.5, 0.025, 0.25], frame_on=False)
            ax.set_xticks([])
            ax.set_xticklabels([])
            ax.set_yticks([])
            ax.set_yticklabels([])
            clb = fig.colorbar(im, ax=ax, fraction=0.4, label='Bootstrap counts corr.')
            clb.ax.tick_params(labelsize=4)

            plt.savefig(figName, dpi=300)

            return

        keys = np.array([key.split('/')[1:] for key in KeysOfStore('results_DCS_35_SRA.h5')])

        for species in ['Mus musculus', 'Homo sapiens', 'SRA']:
            b3, b4 = getDf(keys, 'combo3', 'b', species), getDf(keys, 'combo4', 'b', species)

            b3.to_excel('%s b3.xlsx' % species)
            b4.to_excel('%s b4.xlsx' % species)

            plotDf(b3.corr(), '%s b3.png' % species)
            plotDf(b4.corr(), '%s b4.png' % species)

            r3, r4 = getDf(keys, 'combo3', 'r', species), getDf(keys, 'combo4', 'r', species)
            b3[b3 <= r3.max()] = 0.
            b4[b4 <= r4.max()] = 0.

            plotDf(b3.corr(), '%s b3_cut.png' % species)
            plotDf(b4.corr(), '%s b4_cut.png' % species)
