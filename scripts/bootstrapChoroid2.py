from scRegulation.commonFunctions import *
from scRegulation.analysisPipeline import Analysis

def prepareInputData_human_Choroid_remapped():

    if platform.system() == "Windows":
        fname_all_labels = 'd:/Projects/A_Endothelial/Eye/Voigt_choroid_4567_remapped_bySD_annotation_labels_13700cells.h5'
        fname_all_expr = 'd:/Projects/A_Endothelial/Eye/Voigt_choroid_4567_remapped_bySD_DCS_all.h5'
    else:
        fname_all_labels = '/mnt/research/piermarolab/data/Eye/Voigt_choroid_4567_remapped_bySD_annotation_labels_13700cells.h5'
        fname_all_expr = '/mnt/research/piermarolab/data/Eye/Voigt_choroid_4567_remapped_bySD_DCS_all.h5'

    se_all_labels = pd.read_hdf(fname_all_labels, key='df')

    if False:
        one_index = se_all_labels[(se_all_labels == 'Endothelial') | (se_all_labels == 'Mural')].index
        other_index = se_all_labels.index.difference(one_index)
    else:
        #one_index = pd.concat([se_all_labels[(se_all_labels == 'Endothelial')].sample(n=1600, replace=False), se_all_labels[(se_all_labels == 'Mural')]], axis=0, sort=False).index
        #other_index = se_all_labels[(se_all_labels != 'Endothelial') & (se_all_labels != 'Mural')].index

        one_index = pd.concat([se_all_labels[(se_all_labels == 'Endothelial')].sample(n=1895, replace=False), 
                               se_all_labels[(se_all_labels == 'Macrophage')]], axis=0, sort=False).index
        other_index = se_all_labels[(se_all_labels != 'Endothelial') & 
                                    (se_all_labels != 'Macrophage') &
                                    (se_all_labels != 'Unassigned')].index

    print(one_index.shape[0], dict(zip(*np.unique(se_all_labels[one_index].values, return_counts=True))))
    print(other_index.shape[0], dict(zip(*np.unique(se_all_labels[other_index].values, return_counts=True))))

    df_all = pd.read_hdf(fname_all_expr, key='df')
    print(df_all.shape)

    df_one = df_all[one_index]
    df_one = df_one.loc[df_one.sum(axis=1) > 0.]
    print(df_one)

    df_other = df_all[other_index].reindex(df_one.index)
    print(df_other)

    return df_one, df_other

if __name__ == '__main__':

    if platform.system()=="Windows":
        wdir = 'd:/Projects/A_Endothelial/VS/Endothelial/results/' 
    else:
        wdir = '/mnt/home/domansk6/Projects/Endothelial/results/'

    aHuman = Analysis(workingDir=wdir + 'choroid Voigt remapped/', 
                    otherCaseDir=wdir + 'PanglaoDB_byDCS_mouse/bootstrap/All/',
                    genesOfInterest=receptorsListHugo_2555, 
                    knownRegulators=gEC23, 
                    nCPUs=4 if platform.system()=="Windows" else 20, 
                    panels=['combo3avgs', 'combo4avgs', 'fraction', 'binomial', 'markers', 'top50'], 
                    nBootstrap=100, 
                    perEachOtherCase=False)

    #aHuman.prepareDEG(*prepareInputData_human_Choroid_remapped())

    #aHuman.preparePerBatchCase(exprCutoff=0.01)
    #aHuman.prepareBootstrapExperiments(parallel=True)
    #aHuman.analyzeBootstrapExperiments()
            
    #aHuman.reanalyzeMain()
    #aHuman.analyzeCombinationVariant('Avg combo3avgs')
    #aHuman.analyzeCombinationVariant('Avg combo4avgs')

    aHuman.analyzeAllPeaksOfCombinationVariant('Avg combo3avgs', nG=15, nE=30, fcutoff=0.5, width=50)
    #aHuman.analyzeAllPeaksOfCombinationVariant('Avg combo4avgs', nG=8, nE=30, fcutoff=0.5, width=50)

    #aHuman.scramble(['Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo3/', M=20)
    #aHuman.scramble(['Markers', 'Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo4/', M=20)