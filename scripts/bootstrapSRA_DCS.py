from scRegulation.commonFunctions import *
from scRegulation.analysisPipeline import Analysis

if platform.system() == "Windows":
    wdir = 'd:/Projects/A_Endothelial/VS/Endothelial/results/DCS_SRAs/'
else:
    wdir = '/mnt/research/piermarolab/Sergii/DCS_SRAs/'

if __name__ == '__main__':

    cfile = 'results_DCS_SRAs.h5'

    # Analyze each SRA
    if False:
        cases = pd.read_excel(wdir + 'DCS tissues SRAs.xlsx', sheet_name='Selected DCS', index_col=[0,1,2], header=0)
        cases = cases[cases['cells'] >= 300]
        cases = cases['SRS'].str.replace(' ', '').str.split(',')
        print(cases)

        for (tissue, species, SRA), SRSs in cases.iteritems():
            print('\n', tissue, species, SRA, len(SRSs), flush=True)

            args = dict(genesOfInterest=receptorsListHugo_2555, knownRegulators=gEC23, nCPUs=4 if platform.system()=="Windows" else 10, panels=['combo3avgs', 'combo4avgs', 'fraction', 'binomial', 'markers', 'top50'], nBootstrap=100, perEachOtherCase=False)

            allPanglaoDBmouse = '/mnt/home/domansk6/Projects/Endothelial/results/PanglaoDB_byDCS_mouse/bootstrap/All/'
            allPanglaoDBhuman = '/mnt/home/domansk6/Projects/Endothelial/results/PanglaoDB_byDCS_human/bootstrap/All/'
        
            an = Analysis(**dict(args, workingDir=wdir + 'PanglaoDB_byDCS %s %s %s/' % (tissue, species, SRA), otherCaseDir=allPanglaoDBmouse if species == 'Homo sapiens' else allPanglaoDBhuman))

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
