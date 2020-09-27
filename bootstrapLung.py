from commonFunctions import *
from analysisPipeline import Analysis

def prepareDEGforTissues(dataSaveName, species, tissues, preprocessedDataPath = '/mnt/research/piermarolab/Sergii/Endothelial by PanglaoDB definition/'):

    from PanglaoDBannotation import getAnnotationsSummaryDf
    se_anno = getAnnotationsSummaryDf(MetadataDirName)[['Tissue origin of the sample', 'Species']].droplevel('Cluster index').set_index('Species', append=True)[['Tissue origin of the sample']]
    se_anno.columns = ['tissue']
    se_anno = se_anno.loc[~se_anno.index.duplicated(keep='first')]
    se_anno.index.names = ['SRA', 'SRS', 'species']
    se_anno['batch'] = se_anno.index.get_level_values('SRA') + '_' + se_anno.index.get_level_values('SRS')
    se_anno = se_anno.reorder_levels(['species', 'SRA', 'SRS']).set_index('tissue', append=True)['batch']
    se_anno = se_anno.xs(species, level='species').droplevel(['SRA', 'SRS'])

    if type(tissues) is str:
        tissues = [tissues]

    batches = se_anno[np.isin(se_anno.index.values, tissues)].values

    print('Species:', species)
    print('Tissues:', tissues)

    print('Loading all expression data', flush=True)
    df_EC = pd.read_hdf(preprocessedDataPath + 'PanglaoDB_EC_normed_and_filtered.h5', key=species)
    df_EC = df_EC[df_EC.columns[np.isin(df_EC.columns.get_level_values('batch').values, batches)]]
    df_EC = df_EC.loc[(df_EC != 0).sum(axis=1) > 0]
    print(df_EC, flush=True)

    print('Recording selected batches expression data', flush=True)
    df_EC.to_hdf(an.dataSaveName, key='df', mode='a', complevel=4, complib='zlib')

    print('Loading DE ranks data', flush=True)
    df_ranks = pd.read_hdf(preprocessedDataPath + 'PanglaoDB_ttest_ranks_per_batch %s.h5' % species, key='df')
    df_ranks = df_ranks[df_ranks.columns.intersection(batches)]
    df_ranks = df_ranks.loc[~(df_ranks==-1).all(axis=1)]
    print(df_ranks, flush=True)

    print('Recording selected batches DE ranks data', flush=True)
    df_ranks.to_hdf(an.dataSaveName, key='df_ranks', mode='a', complevel=4, complib='zlib')

    return

if __name__ == '__main__':

    args = dict(genesOfInterest=receptorsListHugo_2555, 
                knownRegulators=gEC23, 
                nCPUs=4 if platform.system()=="Windows" else 8, 
                panels=['combo3avgs', 'combo4avgs', 'fraction', 'binomial', 'markers', 'top50'], 
                nBootstrap=15, 
                perEachOtherCase=False)

    mousePanglaoAllBatchesDir = '/mnt/research/piermarolab/Sergii/Endothelial by PanglaoDB definition/EC all bootstrap 100 w21/Mus musculus/All/'
    humanPanglaoAllBatchesDir = '/mnt/research/piermarolab/Sergii/Endothelial by PanglaoDB definition/EC all bootstrap 100 w21/Homo sapiens/All/'

    an = Analysis(**dict(args, workingDir='results/PanglaoDB_lung_mouse/', otherCaseDir=humanPanglaoAllBatchesDir))

    if not os.path.isfile(an.dataSaveName):
        prepareDEGforTissues(an.dataSaveName, 'Mus musculus', ['Lung', 'Lung mesenchyme', 'Fetal lung', 'Lung endoderm'])
    
    an.preparePerBatchCase(exprCutoff=0.05)
    an.prepareBootstrapExperiments()
    an.analyzeBootstrapExperiments()
    an.reanalyzeMain()

    an.analyzeCombinationVariant('Avg combo3avgs')
    an.analyzeCombinationVariant('Avg combo4avgs')
    an.scramble(['Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo3/', M=10)
    an.scramble(['Markers', 'Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo4/', M=10)  