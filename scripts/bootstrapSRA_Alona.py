from scRegulation.commonFunctions import *
from scRegulation.analysisPipeline import Analysis
from SRAs_Alona_DCS_analysis import plot_analyze_SRAs, umap_SRAs

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

    if True:
        plot_analyze_SRAs(wdir, cfile)
        umap_SRAs(wdir)