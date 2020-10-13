from scRegulation.commonFunctions import *
from scRegulation.analysisPipeline import Analysis

if __name__ == '__main__':

    if platform.system()=="Windows":
        wdir = 'd:/Projects/A_Endothelial/VS/Endothelial/results/Alona_SRAs/'
    else:
        wdir = '/mnt/research/piermarolab/Sergii/Alona_SRAs/'

    # Prepare SRA SRS cells counts
    if False:
            df = pd.concat([pd.read_hdf(wdir + 'KDR.h5', key='Homo sapiens'), pd.read_hdf(wdir + 'KDR.h5', key='Mus musculus')], 
                            keys=['Homo sapiens', 'Mus musculus'], names=['species'], axis=0, sort=False).reset_index()
            df['SRA'] = df['batch'].str.split('_', expand=True)[0]
            df['SRS'] = df['batch'].str.split('_', expand=True)[1]
            df = df.drop(['batch', 'KDR'], axis=1).set_index(['SRA', 'SRS', 'species'])
            df = df.groupby(level=['SRA', 'SRS', 'species']).count()
            df = df.loc[df['cell'] >= 10]
            df.columns = ['cells']
            print(df)

            # New PanglaoDB tissue associations
            se = pd.read_excel(wdir + 'PanglaoDB_tissue_groups GP.xlsx', index_col=1)['New Groups'].str.replace('?', '')
            se[se!=se] = se.index.values[se!=se]
            tissueDict = se.to_dict()

            dft = pd.read_excel('dev/PanglaoDBdata/df_cell_type_annotations.xlsx', index_col=[0,1,2], header=0, sheet_name='tissues')
            dft['tissue'] = dft['tissue'].replace(tissueDict)
            df['tissue'] = dft.loc[df.index]
            df = df.reset_index().set_index(['tissue', 'species'])

            dfc = df.set_index('SRA', append=True)
            dfc['SRS count'] = dfc['SRS'].groupby(['tissue', 'species', 'SRA']).agg('unique').apply(len).reindex(dfc.index)
            print(dfc)

            writer = pd.ExcelWriter(wdir + 'Alona tissues counts S.xlsx')
            dfc.to_excel(writer, 'Non-aggregated', merge_cells=False)

            dfc = dfc.groupby(['tissue', 'species', 'SRA']).agg({'SRS':'unique', 'cells':'sum', 'SRS count':'max'})
            dfc['SRS'] = dfc['SRS'].apply(cleanListString)
            print(dfc)
            dfc.to_excel(writer, 'Aggregated', merge_cells=False)

            df = pd.concat([df['SRA'].groupby(level=['tissue', 'species']).agg(lambda s: np.unique(s).shape[0]),
                            df['SRS'].groupby(level=['tissue', 'species']).agg(lambda s: np.unique(s).shape[0]),
                            df['cells'].groupby(level=['tissue', 'species']).sum()], axis=1, sort=False).unstack('species')
            df.columns.names = ['measure', 'species']
            df = df.reorder_levels(['species', 'measure'], axis=1).sort_index(axis=1)
            print(df)
            df.to_excel(writer, 'Summary', merge_cells=False)

            writer.save()

            exit()

    if True:
        cases = pd.read_excel(wdir + 'Alona tissues counts.xlsx', sheet_name='Selected', index_col=[0,1,2], header=0)
        cases = cases[cases['cells'] >= 100]
        cases = cases['SRS'].str.replace(' ', '').str.split(',')
        print(cases)

        for (tissue, species, SRA), SRSs in cases.iteritems():
            print('\n', tissue, species, SRA, len(SRSs), flush=True)

            args = dict(genesOfInterest=receptorsListHugo_2555, knownRegulators=gEC23, nCPUs=4 if platform.system()=="Windows" else 10, panels=['combo3avgs', 'combo4avgs', 'fraction', 'binomial', 'markers', 'top50'], nBootstrap=100, perEachOtherCase=False)

            allPanglaoDBmouse = '/mnt/research/piermarolab/Sergii/PanglaoDB_byAlona/PanglaoDB_byDCS_mouse/bootstrap/All/'
            allPanglaoDBhuman = '/mnt/research/piermarolab/Sergii/PanglaoDB_byAlona/PanglaoDB_byDCS_human/bootstrap/All/'
        
            an = Analysis(**dict(args, workingDir=wdir + 'Alona %s %s %s/' % (tissue, species, SRA), otherCaseDir=allPanglaoDBmouse if species == 'Homo sapiens' else allPanglaoDBhuman))

            if False:
                if os.path.isfile(an.workingDir + 'results.png'):
                    continue

                if not os.path.isfile(an.dataSaveName):
                    df_EC = pd.concat([pd.read_hdf('/mnt/research/piermarolab/Sergii/PanglaoDB_byAlona/PanglaoDB_Alona_EC.h5', key=SRA + '_' + SRS) for SRS in SRSs], axis=1, sort=False).fillna(0.)
                    df_other = pd.concat([pd.read_hdf('/mnt/research/piermarolab/Sergii/PanglaoDB_byAlona/PanglaoDB_Alona_nonEC.h5', key=SRA + '_' + SRS) for SRS in SRSs], axis=1, sort=False).fillna(0.)

                    if len(SRSs) < 5:
                        nBatches = 10
                        df_EC.columns = pd.MultiIndex.from_arrays([np.random.permutation(np.hstack([np.array([str(i)]*len(v), dtype=str) for i, v in enumerate(np.array_split(df_EC.columns.get_level_values('batch').values, nBatches))])), df_EC.columns.get_level_values('cell')], names=['batch', 'cell'])
                        df_other.columns = pd.MultiIndex.from_arrays([np.random.permutation(np.hstack([np.array([str(i)]*len(v), dtype=str) for i, v in enumerate(np.array_split(df_other.columns.get_level_values('batch').values, nBatches))])), df_other.columns.get_level_values('cell')], names=['batch', 'cell'])

                    print(df_EC, df_other)
                    df_EC = df_EC.T.loc[~df_EC.T.index.duplicated(keep='first')].T
                    df_other = df_other.T.loc[~df_other.T.index.duplicated(keep='first')].T
                    print(df_EC.shape, df_other.shape)

                    an.prepareDEG(df_EC, df_other)

            if False:
                an.preparePerBatchCase(exprCutoff=0.05)
                an.prepareBootstrapExperiments(parallel=True)

            if True:
                an.analyzeBootstrapExperiments()
                an.reanalyzeMain()
                an.analyzeCombinationVariant('Avg combo3avgs')
                an.analyzeAllPeaksOfCombinationVariant('Avg combo3avgs', nG=8, nE=30, fcutoff=0.5, width=50)
                an.analyzeCombinationVariant('Avg combo4avgs')
                an.analyzeAllPeaksOfCombinationVariant('Avg combo4avgs', nG=8, nE=30, fcutoff=0.5, width=50)
                an.scramble(['Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo3/', M=20)  
                an.scramble(['Markers', 'Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo4/', M=20)  

            if False:
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

    if False:
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

        # Plots with UMAP layout
        if False:
            df3 = pd.read_excel('for meeting 10 08 2020/b3.xlsx', index_col=0, header=0)
            print(df3)
            df4 = pd.read_excel('for meeting 10 08 2020/b4.xlsx', index_col=0, header=0)
            print(df4)

            import umap

            coords3 = umap.UMAP(random_state=42).fit_transform(df3.values.T).T
            plt.scatter(coords3[0], coords3[1])
            for i, label in enumerate(df3.columns):
                plt.text(coords3[0][i], coords3[1][i], label)
            plt.show()

            coords4 = umap.UMAP(random_state=42).fit_transform(df4.values.T).T
            plt.scatter(coords4[0], coords4[1])
            for i, label in enumerate(df4.columns):
                plt.text(coords4[0][i], coords4[1][i], label)
            plt.show()
