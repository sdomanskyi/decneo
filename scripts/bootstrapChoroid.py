from commonFunctions import *
from analysisPipeline import Analysis

def prepareInputData_human_Choroid_remapped():

    if platform.system() == "Windows":
        fname_EC_expr = 'd:/Projects/A_Endothelial/Eye/Voigt_choroid_4567_remapped_bySD_DCS_EC.h5'
        fname_all_expr = 'd:/Projects/A_Endothelial/Eye/Voigt_choroid_4567_remapped_bySD_DCS_all.h5'
    else:
        fname_EC_expr = '/mnt/research/piermarolab/data/Eye/Voigt_choroid_4567_remapped_bySD_DCS_EC.h5'
        fname_all_expr = '/mnt/research/piermarolab/data/Eye/Voigt_choroid_4567_remapped_bySD_DCS_all.h5'

    df_EC = pd.read_hdf(fname_EC_expr, key='df')
    df_EC = df_EC.loc[df_EC.sum(axis=1) > 0.]
    print(df_EC)

    df_all = pd.read_hdf(fname_all_expr, key='df')
    df_other = df_all[df_all.columns.difference(df_EC.columns)].reindex(df_EC.index)
    print(df_other)

    return df_EC, df_other

if __name__ == '__main__':

    args = dict(genesOfInterest=receptorsListHugo_2555, 
                knownRegulators=gEC23, 
                nCPUs=8 if platform.system()=="Windows" else 10, 
                panels=['combo3avgs', 'combo4avgs', 'fraction', 'binomial', 'markers', 'top50'], 
                nBootstrap=100, 
                perEachOtherCase=False)

    aHuman = Analysis(**dict(args, workingDir='results/choroid Voigt remapped/', otherCaseDir='results/workflow/PanglaoDB EC all cells w21/Mus musculus/'))

    aHuman.prepareDEG(*prepareInputData_human_Choroid_remapped())
    aHuman.preparePerBatchCase(exprCutoff=0.01)
    aHuman.prepareBootstrapExperiments()

    aHuman.analyzeBootstrapExperiments()
            
    aHuman.reanalyzeMain()
    aHuman.analyzeCombinationVariant('Avg combo3avgs')
    aHuman.analyzeCombinationVariant('Avg combo4avgs')

    aHuman.scramble(['Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo3/', M=20)
    aHuman.scramble(['Markers', 'Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo4/', M=20)

    aHuman.analyzeAllPeaksOfCombinationVariant('Avg combo3avgs', nG=8, nE=30, fcutoff=0.5, width=50)
    aHuman.analyzeAllPeaksOfCombinationVariant('Avg combo4avgs', nG=8, nE=30, fcutoff=0.5, width=50)
