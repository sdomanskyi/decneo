from scRegulation.commonFunctions import *

if __name__ == '__main__':

    wdir = 'd:/Projects/A_Endothelial/VS/Endothelial/results/for meeting 10 15 2020/'

    # Something comparison? What am I looking at?
    if False:
        df = pd.read_excel(wdir + 'sel DCS 3 grouped.xlsx', index_col=[0,1], header=[0,1]).groupby(level=0, axis=1).mean()
        df = df.xs(key=0, level=1, axis=0)
        #df = df.droplevel(1, axis=0)
        df = df[df.sum(axis=1) > 0.]
        df.columns = ['brain', 'lung', 'testis']
        print(df)

        dfAlex = pd.read_excel('d:/Projects/A_Endothelial/VS/Endothelial/dev/PL_Norm_Peak_diff SD marked.xlsx', index_col=0, header=0)
        dfAlex.columns = ['brain', 'lung', 'testis']
        dfAlex = dfAlex.loc[df.index]
        print(dfAlex)

        print(df.corrwith(dfAlex))

    # Alona vs DCS all
    if False:
        def getM(df, s):
            df = pd.concat([df[df.columns[0:2]].set_index('3m').dropna(), 
                            df[df.columns[2:4]].set_index('4m').dropna(), 
                            df[df.columns[4:6]].set_index('3h').dropna(), 
                            df[df.columns[6:8]].set_index('4h').dropna()], 
                            axis=1, sort=False, keys=['3m' + s, '4m' + s, '3h' + s, '4h' + s]).fillna(0.).droplevel(1, axis=1)

            return df

        dfD = getM(pd.read_excel(wdir + 'summary all PanglaoDB.xlsx', sheet_name='DCS', inde_col=None, header=0), '_DCS')
        dfA = getM(pd.read_excel(wdir + 'summary all PanglaoDB.xlsx', sheet_name='Alona', inde_col=None, header=0), '_Alona')
        df = pd.concat([dfA, dfD], axis=1, sort=False).fillna(0.)
        print(df)

        df.to_excel(wdir + 'summary_aligned.xlsx')

    # Part main
    if False:
        def getDfC(tool, combo, suff):
            df = pd.read_excel(wdir + '%s_SRAs/' % tool + 'SRA %s%s.xlsx' % (suff, combo), index_col=0, header=0)
        
            df.loc['combo'] = [combo for c in df.columns]
            df.loc['tool'] = [tool for c in df.columns]
            df.loc['species'] = ['Homo sapiens' if 'Homo sapiens' in c else 'Mus musculus' for c in df.columns]
            df.loc['tissue'] = [c.split(' Homo sapiens ')[0] if 'Homo sapiens' in c else c.split(' Mus musculus ')[0] for c in df.columns]
            df = df.T.set_index('combo', append=True).T
            df = df.T.set_index('tool', append=True).T
            df = df.T.set_index('species', append=True).T
            df = df.T.set_index('tissue', append=True).T
            df = df.astype(float)
            df.columns.names = ['id', 'combo', 'tool', 'species', 'tissue']
            df = df.groupby(level=['combo', 'tool', 'species', 'tissue'], axis=1).mean()

            return df

        dfr = pd.concat([getDfC('Alona', 3, 'r'), 
                        getDfC('DCS', 3, 'r'),
                        getDfC('Alona', 4, 'r'), 
                        getDfC('DCS', 4, 'r')], axis=1).fillna(0.).max(axis=0).to_frame().T
        df = pd.concat([getDfC('Alona', 3, 'b'), 
                        getDfC('DCS', 3, 'b'),
                        getDfC('Alona', 4, 'b'), 
                        getDfC('DCS', 4, 'b')], axis=1).fillna(0.)
        df[df.values <= dfr.values] = 0.
        df = df.loc[df.sum(axis=1) > 0.]

        dfa = pd.read_excel(wdir + 'summary all PanglaoDB aligned.xlsx', sheet_name='Sheet1', index_col=0, header=0)
        dfa.columns = pd.MultiIndex.from_tuples([(3, 'Alona', 'Mus musculus', 'All'),
                                                 (4, 'Alona', 'Mus musculus', 'All'),
                                                 (3, 'Alona', 'Homo sapiens', 'All'),
                                                 (4, 'Alona', 'Homo sapiens', 'All'),
                                                 (3, 'DCS', 'Mus musculus', 'All'),
                                                 (4, 'DCS', 'Mus musculus', 'All'),
                                                 (3, 'DCS', 'Homo sapiens', 'All'),
                                                 (4, 'DCS', 'Homo sapiens', 'All')], names=['combo', 'tool', 'species', 'tissue'])
        dfa[dfa <= np.array([0.238, 0.239, 0.378, 0.279, 0.244, 0.240, 0.292, 0.254])] = 0.
        dfa = dfa.loc[dfa.sum(axis=1) > 0.]

        dfc3 = pd.read_excel(wdir + 'choroid_remapped.xlsx', sheet_name='3', index_col=0, header=[0,1,2,3])
        dfc3[dfc3 <= 0.2579] = 0.
        dfc3 = dfc3.loc[dfc3.sum(axis=1) > 0.]
        dfc4 = pd.read_excel(wdir + 'choroid_remapped.xlsx', sheet_name='4', index_col=0, header=[0,1,2,3])
        dfc4[dfc4 <= 0.25799] = 0.
        dfc4 = dfc4.loc[dfc4.sum(axis=1) > 0.]

        df = pd.concat([dfa, df, dfc3, dfc4], axis=1, sort=False).fillna(0.).astype(float)
        df.T.to_excel(wdir + 'counts_cut.xlsx', merge_cells=False)

        df = df.T

        if False:
            writer = pd.ExcelWriter(wdir + 'z.xlsx')
            for combo, tool in [(3, 'Alona'), (3, 'DCS'), (4, 'Alona'), (4, 'DCS')]:
                df_temp = df.xs(key=combo, level='combo').xs(key=tool, level='tool')
                df_temp = df_temp[df_temp.columns[df_temp.sum(axis=0) > 0.]]
                df_temp[:] = (df_temp.values - np.mean(df_temp.values, axis=0)) / np.std(df_temp.values, axis=0)
                df_temp = df_temp[df_temp.columns[np.abs(df_temp).sum(axis=0) > 0.]]

                df_temp = df_temp.iloc[scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(df_temp, 'ward'), no_plot=True, get_leaves=True)['leaves']]
                df_temp = df_temp.T.iloc[scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(df_temp.T, 'ward'), no_plot=True, get_leaves=True)['leaves']].T

                df_temp.to_excel(writer, sheet_name='combo%s %s' % (combo, tool), merge_cells=False)
                print(df_temp)

                fig, ax = plt.subplots(figsize=(10, 4))
                ax.imshow(df_temp.values, cmap='Greens', interpolation='None', aspect='auto')

                ax.set_xticks(range(df_temp.shape[1]))
                ax.set_xticklabels(df_temp.columns, rotation=90, fontsize=3)
                ax.set_xlim([-0.5, df_temp.shape[1] - 0.5])

                ax.set_yticks(range(df_temp.shape[0]))
                ax.set_yticklabels(['%s (%s)' % (col[1], col[0]) for col in df_temp.index], rotation=0, fontsize=5)
                ax.set_ylim([-0.5, df_temp.shape[0] - 0.5])

                ax.set_title('combo%s %s' % (combo, tool), fontdict={'color': 'b', 'size':'8'})
                fig.tight_layout()
                fig.savefig(wdir + 'combo%s %s.png' % (combo, tool), dpi=600)
                fig.savefig(wdir + 'combo%s %s.pdf' % (combo, tool))

            writer.save()

        if False:
            writer = pd.ExcelWriter(wdir + 'figure data bootstrap counts.xlsx')

            for combo, tool in [(3, 'Alona'), (3, 'DCS'), (4, 'Alona'), (4, 'DCS')]:
                df_temp = df.xs(key=combo, level='combo').xs(key=tool, level='tool')
                df_temp = df_temp[df_temp.columns[df_temp.sum(axis=0) > 0.]]

                #df_temp[:] = (df_temp.values - np.mean(df_temp.values, axis=0)) / np.std(df_temp.values, axis=0)
                #df_temp = df_temp[df_temp.columns[np.abs(df_temp).sum(axis=0) > 0.]]

                df_temp = df_temp.iloc[scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(df_temp, 'ward'), no_plot=True, get_leaves=True)['leaves']]
                df_temp = df_temp.T.iloc[scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(df_temp.T, 'ward'), no_plot=True, get_leaves=True)['leaves']].T

                df_temp.to_excel(writer, sheet_name='combo%s %s' % (combo, tool), merge_cells=False)
                print(df_temp)

                fig, ax = plt.subplots(figsize=(10, 4))
                ax.imshow(df_temp.values, cmap='Greens', interpolation='None', aspect='auto')

                ax.set_xticks(range(df_temp.shape[1]))
                ax.set_xticklabels(df_temp.columns, rotation=90, fontsize=2)
                ax.set_xlim([-0.5, df_temp.shape[1] - 0.5])

                ax.set_yticks(range(df_temp.shape[0]))
                ax.set_yticklabels(['%s (%s)' % (col[1], col[0]) for col in df_temp.index], rotation=0, fontsize=5)
                ax.set_ylim([-0.5, df_temp.shape[0] - 0.5])

                ax.set_title('combo%s %s' % (combo, tool), fontdict={'color': 'b', 'size':'8'})
                fig.tight_layout()
                fig.savefig(wdir + 'combo%s %s bootstrap counts.png' % (combo, tool), dpi=600)
                fig.savefig(wdir + 'combo%s %s bootstrap counts.pdf' % (combo, tool))

            writer.save()

    # Select tissue-specific genes DCS / Alona
    if False:
        g0 = []
        g1 = []
        for tool in ['DCS', 'Alona']:
            df = pd.read_excel(wdir + 'sel %s 3.xlsx' % tool, index_col=[0,1], header=[0,1])
            #print(df)

            se0s = []
            se1s = []
            for col in df.columns:
                se = df[col].sort_values(ascending=False)
                se0 = se.xs(0, level=1)[:100]
                se0 = se0[se0 >= 0.25]
                se1 = se.xs(1, level=1)
                se1 = se1[se1 < 0.25]

                se0 = pd.Series(se0.index)
                se1 = pd.Series(se1.index)

                se0.name = col
                se1.name = col

                se0s.append(se0)
                se1s.append(se1)
        
            se0s = pd.concat(se0s, axis=1, sort=False)
            se1s = pd.concat(se1s, axis=1, sort=False)

            #print(se0s)
            #print(se1s)
            #print()

            g0.append(se0s)
            g1.append(se1s)

            writer = pd.ExcelWriter(wdir + '3 %s top-bottom bootstrap count 10 20 2020.xlsx' % tool)
            se0s.T.to_excel(writer, sheet_name='top', merge_cells=False)
            se1s.T.to_excel(writer, sheet_name='bottom', merge_cells=False)
            writer.save()

        def JaccardOfDataFrames(g0):

            df0 = pd.DataFrame(index=g0[0].columns, columns=g0[1].columns)
            for col0 in g0[0].columns:
                for col1 in g0[1].columns:
                    a, b = set(g0[0][col0].dropna().values[:10]), set(g0[1][col1].dropna().values[:10])
                    u = a.union(b)

                    if len(u) > 0:
                        J = len(a.intersection(b)) / len(u)

                        df0.loc[col0, col1] = J

            print(df0)

            return df0

        writer = pd.ExcelWriter(wdir + 'DCS vs Alona combo3 tissue-specific.xlsx')
        JaccardOfDataFrames(g0).to_excel(writer, sheet_name='Up')
        JaccardOfDataFrames(g1).to_excel(writer, sheet_name='Down')
        writer.save()

    # Part sub-peaks
    if False:
        def getDfC(tool, combo, suff):

            df = pd.read_excel(wdir + '%s_SRAs/' % tool + 'SRA %s%s.xlsx' % (suff, combo), index_col=0, header=0)
            #df['frequency'] = [str(i) + ':' + f for i, f in enumerate(df['frequency'].astype(str))]
            df = df.set_index('frequency', append=True)
            df['combo'] = [combo for c in df.index]
            df['tool'] = [tool for c in df.index]
            df = df.set_index('combo', append=True)
            df = df.set_index('tool', append=True)
            df.index.names = ['id', 'frequency', 'combo', 'tool']
            df['peaks'] = df['peaks'].str.split(', ')

            return df

        df = pd.concat([getDfC('Alona', 3, 's'), 
                        getDfC('DCS', 3, 's'),
                        getDfC('Alona', 4, 's'), 
                        getDfC('DCS', 4, 's')], axis=0)
        print(df)

        writer = pd.ExcelWriter(wdir + 'peaks.xlsx')
        for combo, tool in [(3, 'Alona'), (3, 'DCS'), (4, 'Alona'), (4, 'DCS')]:
            df_temp = df.xs(key=combo, level='combo').xs(key=tool, level='tool')
            df_temp['peaks'] = df_temp['peaks'].apply(np.array)
            print(df_temp)

            # Add all
            if True:
                dfh = pd.read_excel(wdir + '/all peaks/%s combo%s peaks Homo sapiens.xlsx' % (tool, combo), index_col=0, header=0)
                dfh = dfh['genes'].str.split(', ').apply(np.array)
                dfh = pd.concat([dfh], keys=['All Homo sapiens'], names=['id'], sort=False)
                dfh.name = 'peaks'

                dfm = pd.read_excel(wdir + '/all peaks/%s combo%s peaks Mus musculus.xlsx' % (tool, combo), index_col=0, header=0)
                dfm = dfm['genes'].str.split(', ').apply(np.array)
                dfm = pd.concat([dfm], keys=['All Mus musculus'], names=['id'], sort=False)
                dfm.name = 'peaks'

                df_temp = pd.concat([dfm.to_frame(), dfh.to_frame(), df_temp], axis=0, sort=False)

            # Add choroid
            if True:
                if tool == 'DCS':
                    dfch = pd.read_excel(wdir + 'choroid All peaks Avg combo%savgs. nG8-nE30.xlsx' % combo, index_col=0, header=0)
                    dfch = dfch['genes'].str.split(', ').apply(np.array)
                    dfch = pd.concat([dfch], keys=['Choroid Homo sapiens'], names=['id'], sort=False)
                    dfch.name = 'peaks'

                    df_temp = pd.concat([df_temp, dfch.to_frame()], axis=0, sort=False)
            
            print(df_temp)

            genes = np.unique(np.hstack(df_temp['peaks'].values))
            ids = np.unique(df_temp.index.get_level_values('id'))
            print(df_temp.shape, len(genes), len(ids))

            dfm = pd.DataFrame(index=ids, columns=genes, data=0.)
            for id in ids:
                dfid = df_temp.xs(id, level='id')
                for i, v in enumerate(dfid['peaks'].values):
                    dfm.loc[id, v] = dfid.index[i]

            df_temp = dfm
            df_temp = df_temp[df_temp.columns[df_temp.sum(axis=0) > 0.2]]

            df_temp = df_temp.T.iloc[scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(df_temp.T, 'ward'), no_plot=True, get_leaves=True)['leaves']].T
             
            df_temp.to_excel(writer, sheet_name='combo%s %s' % (combo, tool), merge_cells=False)

            # Make plot
            if True:
                fig, ax = plt.subplots(figsize=(10, 4))

                cmap = cm.jet
                cmap.set_bad('white')

                ax.imshow(np.ma.array(df_temp.values, mask=(df_temp.values==0.)), cmap=cmap, interpolation='None', aspect='auto')

                ax.set_xticks(range(df_temp.shape[1]))
                ax.set_xticklabels(df_temp.columns, rotation=90, fontsize=2.5)
                ax.set_xlim([-0.5, df_temp.shape[1] - 0.5])

                ax.set_yticks(range(df_temp.shape[0]))
                ax.set_yticklabels(df_temp.index, rotation=0, fontsize=5)
                ax.set_ylim([-0.5, df_temp.shape[0] - 0.5])

                ax.set_title('combo%s %s' % (combo, tool), fontdict={'color': 'b', 'size':'8'})
                fig.tight_layout()
                fig.savefig(wdir + 'peaks combo%s %s.png' % (combo, tool), dpi=600)
                fig.savefig(wdir + 'peaks combo%s %s.pdf' % (combo, tool))

        writer.save()

    # Part sub-peaks grouped DCS 3
    if True:
        df = pd.read_excel(wdir + 'sel DCS 3.xlsx', index_col=[0,1], header=[0,1])
        df.columns.names = ['species', 'tissue']
        df.index.names = ['gene', 'isGeneric']
        df = df.xs(0, level='isGeneric')
        df = df[df.columns[df.columns.get_level_values('tissue') != 'All']]
        print(df)

        ses1 = []
        for col in df.columns:
            se = df[col][df[col] > 0.1].sort_values(ascending=False)[:15]
            df.loc[:, col] = 0.
            df.loc[se.index, col] = se.values
            se.name = col
            se[:] = se.index.values
            se.index = range(len(se))
            ses1.append(se)

        ses1 = pd.concat(ses1, axis=1, sort=False)
        ses1.to_excel(wdir + 'DCS 3 peaks REGULAR.xlsx', merge_cells=False)
        print(ses1)


        df = pd.read_excel(wdir + 'peaks DCS 3.xlsx', sheet_name='Sheet1', index_col=[0,1], header=[0,1,2])
        df.columns.names = ['tissue', 'species', 'id']
        df.index.names = ['gene', 'isGeneric']
        df = df.groupby(by=['species', 'tissue'], axis=1).mean()
        df = df.xs(0, level='isGeneric')
        df = df[df.columns[df.columns.get_level_values('tissue') != 'All']]
        print(df.shape)

        ses2 = []
        for col in df.columns:
            se = df[col][df[col] > 0.1].sort_values(ascending=False)[:15]
            df.loc[:, col] = 0.
            df.loc[se.index, col] = se.values
            se.name = col
            se[:] = se.index.values
            se.index = range(len(se))
            ses2.append(se)

        ses2 = pd.concat(ses2, axis=1, sort=False)
        ses2.to_excel(wdir + 'DCS 3 peaks chopped lists.xlsx', merge_cells=False)
        print(ses2)


        def JaccardOfDataFrames(g0):

            df0 = pd.DataFrame(index=g0[0].columns, columns=g0[1].columns)
            for col0 in g0[0].columns:
                for col1 in g0[1].columns:
                    a, b = set(g0[0][col0].dropna().values[:10]), set(g0[1][col1].dropna().values[:10])
                    u = a.union(b)
                    if len(u) > 0:
                        df0.loc[col0, col1] = len(a.intersection(b)) / len(u)
            print(df0)

            return df0

        dfj = JaccardOfDataFrames((ses1, ses2))
        dfj = dfj.sort_index(axis=0).sort_index(axis=1)
        dfj.to_excel(wdir + 'DCS 3 regular vs peaks.xlsx', merge_cells=True)


        exit()

        df = df.loc[df.sum(axis=1) > 0.]
        print(df)

        df_temp = df.T

        df_temp = df_temp.iloc[scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(df_temp, 'ward'), no_plot=True, get_leaves=True)['leaves']]
        df_temp = df_temp.T.iloc[scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(df_temp.T, 'ward'), no_plot=True, get_leaves=True)['leaves']].T

        #df_temp.to_excel(wdir + 'DCS 3 peaks chopped.xlsx')



        
        # Make plot
        if False:
            fig, ax = plt.subplots(figsize=(12, 4))

            cmap = cm.jet
            cmap.set_bad('white')

            ax.imshow(np.ma.array(df_temp.values, mask=(df_temp.values==0.)), cmap=cmap, interpolation='None', aspect='auto')

            ax.set_xticks(range(df_temp.shape[1]))
            ax.set_xticklabels(df_temp.columns, rotation=90, fontsize=5)
            ax.set_xlim([-0.5, df_temp.shape[1] - 0.5])

            ax.set_yticks(range(df_temp.shape[0]))
            ax.set_yticklabels(df_temp.index, rotation=0, fontsize=5)
            ax.set_ylim([-0.5, df_temp.shape[0] - 0.5])

            fig.tight_layout()
            fig.savefig(wdir + 'DCS 3 peaks chopped.png', dpi=600)
