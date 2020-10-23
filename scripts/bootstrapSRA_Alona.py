from scRegulation.commonFunctions import *
from scRegulation.analysisPipeline import Analysis

if platform.system() == "Windows":
    wdir = 'd:/Projects/A_Endothelial/VS/Endothelial/results/Alona_SRAs/'
else:
    wdir = '/mnt/research/piermarolab/Sergii/Alona_SRAs/'

if __name__ == '__main__':

    cfile = 'results_Alona_SRAs.h5'

    # Select SRAs to analyze
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

    # Analyze each SRA
    if False:
        cases = pd.read_excel(wdir + 'Alona tissues counts.xlsx', sheet_name='Selected', index_col=[0,1,2], header=0)
        cases = cases[cases['cells'] >= 300]
        cases = cases['SRS'].str.replace(' ', '').str.split(',')
        print(cases)

        ##cases = cases.iloc[:10]
        ##cases = cases.iloc[10:20]
        ##cases = cases.iloc[20:30]
        ##cases = cases.iloc[30:40]
        #cases = cases.iloc[40:]

        for (tissue, species, SRA), SRSs in cases.iteritems():
            print('\n', tissue, species, SRA, len(SRSs), flush=True)

            args = dict(genesOfInterest=receptorsListHugo_2555, knownRegulators=gEC23, nCPUs=4 if platform.system()=="Windows" else 10, panels=['combo3avgs', 'combo4avgs', 'fraction', 'binomial', 'markers', 'top50'], nBootstrap=100, perEachOtherCase=False)

            allPanglaoDBmouse = '/mnt/research/piermarolab/Sergii/PanglaoDB_byAlona/PanglaoDB_byDCS_mouse/bootstrap/All/'
            allPanglaoDBhuman = '/mnt/research/piermarolab/Sergii/PanglaoDB_byAlona/PanglaoDB_byDCS_human/bootstrap/All/'
        
            an = Analysis(**dict(args, workingDir=wdir + 'Alona %s %s %s/' % (tissue, species, SRA), otherCaseDir=allPanglaoDBmouse if species == 'Homo sapiens' else allPanglaoDBhuman))

            # Prepare expression and DEG data
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

            # Run analysis pipeline
            if False:
                an.preparePerBatchCase(exprCutoff=0.05)
                an.prepareBootstrapExperiments(parallel=True)
                an.analyzeBootstrapExperiments()
                an.reanalyzeMain()
                an.analyzeCombinationVariant('Avg combo3avgs')
                an.analyzeAllPeaksOfCombinationVariant('Avg combo3avgs', nG=8, nE=30, fcutoff=0.5, width=50)
                an.analyzeCombinationVariant('Avg combo4avgs')
                an.analyzeAllPeaksOfCombinationVariant('Avg combo4avgs', nG=8, nE=30, fcutoff=0.5, width=50)
                an.scramble(['Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo3/', M=20)  
                an.scramble(['Markers', 'Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo4/', M=20)  

            # Collect results
            if True:
                name = '%s %s %s' % (tissue, species, SRA)

                # Collect bootstrap counts for combo3 and combo4
                b3 = pd.read_excel(an.workingDir + 'Avg combo3avgs_variant.xlsx', index_col=0)['Bootstrap']
                b4 = pd.read_excel(an.workingDir + 'Avg combo4avgs_variant.xlsx', index_col=0)['Bootstrap']
                b3.name = 'combo3'
                b4.name = 'combo4'
                b3.to_hdf(wdir + cfile, key='combo3/b/%s' % name, mode='a', complevel=4, complib='zlib')
                b4.to_hdf(wdir + cfile, key='combo4/b/%s' % name, mode='a', complevel=4, complib='zlib')

                # Collect max random values for combo3 and combo4
                r3 = pd.read_excel(an.workingDir + 'random/combo3/se_distribution.xlsx', index_col=0)[0]
                r4 = pd.read_excel(an.workingDir + 'random/combo4/se_distribution.xlsx', index_col=0)[0]
                r3.name = 'combo3'
                r4.name = 'combo4'
                r3.to_hdf(wdir + cfile, key='combo3/r/%s' % name, mode='a', complevel=4, complib='zlib')
                r4.to_hdf(wdir + cfile, key='combo4/r/%s' % name, mode='a', complevel=4, complib='zlib')

                # Collect sub-peaks lists for combo3 and combo4
                s3 = pd.read_excel(an.workingDir + 'All peaks Avg combo3avgs. nG8-nE30.xlsx', index_col=0)['genes'].str.replace(' ', '').str.split(',')
                s4 = pd.read_excel(an.workingDir + 'All peaks Avg combo4avgs. nG8-nE30.xlsx', index_col=0)['genes'].str.replace(' ', '').str.split(',')
                s3.name = 'combo3'
                s3.to_hdf(wdir + cfile, key='combo3/s/%s' % name, mode='a', complevel=4, complib='zlib')
                s4.to_hdf(wdir + cfile, key='combo4/s/%s' % name, mode='a', complevel=4, complib='zlib')



    # Analyze all the results
    if True:
        def getDf(keys, combo, parameter, species, peaks=False):

            dfs = []
            for name in np.unique(keys.T[2]):
                if species in name:
                    b = pd.read_hdf(wdir + cfile, key='/%s/%s/%s' % (combo, parameter, name))

                    if peaks:
                        b.name = 'peaks'
                        b = pd.concat([b.to_frame()], axis=0, keys=[name])
                    else:
                        b.name = name

                    dfs.append(b)

            if peaks:
                dfs = pd.concat(dfs, axis=0, sort=False).fillna(0.)
                dfs['peaks'] = dfs['peaks'].apply(cleanListString)
            else:
                dfs = pd.concat(dfs, axis=1, sort=False).fillna(0.)
            print(dfs.shape)

            return dfs

        def plotDf(df, figName, n_clusters = 10):
            fig = plt.figure(figsize=(8,10))

            ax = fig.add_axes([0.35, 0.7, 0.55, 0.15], frame_on=False)

            data = np.nan_to_num(df.values, nan=df.values.max())
            Z = hierarchy.linkage(data, method='ward', optimal_ordering=True)
            origLineWidth = matplotlib.rcParams['lines.linewidth']
            matplotlib.rcParams['lines.linewidth'] = 0.5
            # color_threshold = (Z[-n_clusters,2] + Z[-n_clusters+1,2]) / 2
            D = hierarchy.dendrogram(Z, ax=ax, color_threshold = 0, above_threshold_color='k', orientation='top')
            matplotlib.rcParams['lines.linewidth'] = origLineWidth
            ax.set_xticklabels([])
            ax.set_xticks([])
            ax.set_yticklabels([])
            ax.set_yticks([])

            tissues = [c.split(' Homo sapiens ')[0] if 'Homo sapiens' in c else c.split(' Mus musculus ')[0] for c in df.columns]
            dcolors = {c[1]:c[0] for c in list(enumerate(np.unique(tissues)))}
            colors = {k: cm.jet(v/len(dcolors)) for k, v in dcolors.items()}
            colors = [colors[t] for t in tissues]
            colors = np.array(colors)[D['leaves']]

            #cluster_labels = scipy.cluster.hierarchy.fcluster(Z, t=n_clusters, criterion='maxclust')[D['leaves']] - 1
            cluster_labels = np.array(tissues)[D['leaves']]

            #for label in np.unique(cluster_labels):
            #    pos = 5. + 10. * np.where(cluster_labels==label)[0].mean()
            #    ax.text(pos, 0.3, str(label), color='blue', fontsize=14, va='center', ha='center').set_path_effects([path_effects.Stroke(linewidth=1.5, foreground='grey'), path_effects.Normal()])

            #plt.title(figName)


            ax = fig.add_axes([0.35, 0.3, 0.55, 0.4], frame_on=False)
            cmap = plt.cm.hot
            cmap.set_bad('grey')
            sdata = np.ma.array(df.values[:,D['leaves']][D['leaves'],:], mask=np.isnan(df.values))
            im = ax.imshow(sdata, cmap=cmap, aspect='auto', interpolation='None', extent=(-0.5, df.shape[0] - 0.5, df.shape[1] - 0.5, -0.5))
            ax.set_xticks(range(len(df.columns)))
            ax.set_yticks(range(len(df.columns)))
            ax.set_xticklabels(df.columns.values[D['leaves']])
            ax.set_yticklabels(df.columns.values[D['leaves']])
            ax.tick_params(axis='x', labelsize=8, width=0.25, length=1, rotation=90)
            ax.tick_params(axis='y', labelsize=8, width=0.25, length=1, rotation=0)

            for i, tick in enumerate(plt.gca().get_xticklabels()):
                tick.set_color(colors[i])

            for i, tick in enumerate(plt.gca().get_yticklabels()):
                tick.set_color(colors[i])

            ax = fig.add_axes([0.9, 0.5, 0.025, 0.25], frame_on=False)
            ax.set_xticks([])
            ax.set_xticklabels([])
            ax.set_yticks([])
            ax.set_yticklabels([])
            clb = fig.colorbar(im, ax=ax, fraction=0.4, label='Bootstrap counts corr.')
            clb.ax.tick_params(labelsize=4)

            plt.savefig(figName, dpi=300)
            plt.clf()
            plt.close()

            return df.values[:,D['leaves']][D['leaves'],:], n_clusters, cluster_labels

        keys = np.array([key.split('/')[1:] for key in KeysOfStore(wdir + cfile)])
        for species in ['Mus musculus', 'Homo sapiens', 'SRA'][-1:]:
            if False:
                s3, s4 = getDf(keys, 'combo3', 's', species, True), getDf(keys, 'combo4', 's', species, True)
                s3.to_excel(wdir + '%s s3.xlsx' % species, merge_cells=False)
                s4.to_excel(wdir + '%s s4.xlsx' % species, merge_cells=False)

            if False:
                r3, r4 = getDf(keys, 'combo3', 'r', species), getDf(keys, 'combo4', 'r', species)
                r3.to_excel(wdir + '%s r3.xlsx' % species)
                r4.to_excel(wdir + '%s r4.xlsx' % species)

            if True:
                b3, b4 = getDf(keys, 'combo3', 'b', species), getDf(keys, 'combo4', 'b', species)
                b3.to_excel(wdir + '%s b3.xlsx' % species)
                b4.to_excel(wdir + '%s b4.xlsx' % species)
                c3 = plotDf(b3.corr(), wdir + '%s b3.png' % species, n_clusters=6)
                c4 = plotDf(b4.corr(), wdir + '%s b4.png' % species, n_clusters=6)

                silhouette(*c3)

                exit()

    # Plots with UMAP layout
    if False:
        df3 = pd.read_excel(wdir + 'SRA b3.xlsx', index_col=0, header=0)
        df4 = pd.read_excel(wdir + 'SRA b4.xlsx', index_col=0, header=0)

        tissues = [c.split(' Homo sapiens ')[0] if 'Homo sapiens' in c else c.split(' Mus musculus ')[0] for c in df3.columns]
        dcolors = {c[1]:c[0] for c in list(enumerate(np.unique(tissues)))}
        colors = {k: cm.jet(v/len(dcolors)) for k, v in dcolors.items()}
        colors = [colors[t] for t in tissues]

        def pplot(coords, columns, saveName):
            
            plt.scatter(coords[0], coords[1])

            for i, label in enumerate(columns):
                plt.text(coords[0][i], coords[1][i], label, color=colors[i], fontsize=4).set_path_effects([path_effects.Stroke(linewidth=0.1, foreground='grey'), path_effects.Normal()])

            plt.axis('off')
            plt.savefig(saveName, dpi=300)
            plt.clf()

            return

        import umap

        pplot(umap.UMAP(random_state=42).fit_transform(df3.values.T).T, df3.columns, wdir + 'UMAP_b3.png')
        pplot(umap.UMAP(random_state=42).fit_transform(df4.values.T).T, df4.columns, wdir + 'UMAP_b4.png')
