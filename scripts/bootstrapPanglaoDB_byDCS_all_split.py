from decneo.commonFunctions import *
from decneo.analysisPipeline import Analysis

if __name__ == '__main__':
  
    if platform.system()=="Windows":
        wdir = 'd:/Projects/A_Endothelial/VS/Endothelial/results/' 
    else:
        wdir = '/mnt/research/piermarolab/Sergii/results/'

    # Split batches into 2 halfs (randomly)
    if False:
        df_measure = pd.read_hdf(wdir + 'PanglaoDB_byDCS_mouse_correlation/byBatches/metricsFile.h5')
        print(df_measure)

        np.random.seed(0)
        if True:
            # Split by SRAs
            SRAs = df_measure.columns.str.split('_', expand=True).get_level_values(0)
            uSRAs = np.unique(SRAs)
            selSRAs = np.random.choice(uSRAs, size=int(len(uSRAs)/2), replace=False)
            print(selSRAs)

            batches = df_measure.columns[np.isin(SRAs, selSRAs)]
        else:
            # Split by SRSs
            batches = np.random.choice(df_measure.columns, size=int(len(df_measure.columns)/2), replace=False)

        print(len(batches))

        print('Recording split 1')
        df_measure[df_measure.columns[np.isin(df_measure.columns, batches)]].to_hdf(wdir + 'PanglaoDB_byDCS_mouse_correlation_splitSRA1/byBatches/metricsFile.h5', key='correlation', **phdf)

        print('Recording split 2')
        df_measure[df_measure.columns[~np.isin(df_measure.columns, batches)]].to_hdf(wdir + 'PanglaoDB_byDCS_mouse_correlation_splitSRA2/byBatches/metricsFile.h5', key='correlation', **phdf)

    # Process both halfs
    if False:
        def doSplitPart(ind):

            an = Analysis(workingDir=wdir + 'PanglaoDB_byDCS_mouse_correlation_splitSRA%s/' % ind, 
                          otherCaseDir=wdir + 'PanglaoDB_byDCS_human_correlation/', 
                          genesOfInterest=receptorsListHugo_2555, knownRegulators=gEC22, nCPUs=4 if platform.system()=="Windows" else 15,
                          perEachOtherCase=True, majorMetric='correlation', dendrogramMetric='euclidean', dendrogramLinkageMethod='ward')

            an.prepareBootstrapExperiments(parallel=False)
            an.analyzeBootstrapExperiments()
            an.reanalyzeMain()
            an.analyzeCombinationVariant('Avg combo3avgs')
            an.analyzeCombinationVariant('Avg combo4avgs')
            an.generateAnalysisReport()

            return

        doSplitPart(1) # ~4 hrs, 100 Gb
        doSplitPart(2) # ~4 hrs, 100 Gb

    # Compare the two halfs
    if True:
        subDir = 'results majorMetric=correlation, dendrogramMetric=euclidean, linkageMethod=ward/'
        dirs = [wdir + 'PanglaoDB_byDCS_mouse_correlation/' + subDir,
                wdir + 'PanglaoDB_byDCS_mouse_correlation_splitSRA1/' + subDir,
                wdir + 'PanglaoDB_byDCS_mouse_correlation_splitSRA2/' + subDir]
        names = ['all', 'half1', 'half2']

        df = []
        for dir, name in zip(dirs, names):
            se = pd.read_excel(dir + 'Avg combo3avgs_variant.xlsx', index_col=0, header=0)['Bootstrap']
            se.name = name
            df.append(se)

        df = pd.concat(df, axis=1, sort=False)
        df = df.sort_values(by='all', ascending=False)
        print(df)
        df.to_excel(wdir + 'all_half1_half2_SRA.xlsx')
