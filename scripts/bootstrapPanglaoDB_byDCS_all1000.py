from decneo.commonFunctions import *
from decneo.analysisPipeline import Analysis, process

def process_loc(df1main, df1other, df2main, df2other, dir1, dir2, genesOfInterest = None, knownRegulators = None, nCPUs = 4,
            panels = ['fraction', 'binomial', 'top50', 'markers', 'combo3avgs', 'combo4avgs'], parallelBootstrap = False,
            exprCutoff1 = 0.05, exprCutoff2 = 0.05, perEachOtherCase = True, doScramble = False, part1 = True, part2 = True, bootstrapDir = 'bootstrap/', **kwargs):

    if genesOfInterest is None:
        genesOfInterest = receptorsListHugo_2555
    
    if knownRegulators is None:
        knownRegulators = gEC22

    an1 = Analysis(workingDir=dir1, otherCaseDir=dir2, genesOfInterest=genesOfInterest,
                        knownRegulators=knownRegulators, panels=panels, nCPUs=nCPUs, perEachOtherCase=perEachOtherCase, **kwargs)

    an1.bootstrapDir = an1.workingDir + bootstrapDir

    if perEachOtherCase:
        an2 = Analysis(workingDir=dir2, otherCaseDir=dir1, genesOfInterest=genesOfInterest,
                            knownRegulators=knownRegulators, panels=panels, nCPUs=nCPUs, perEachOtherCase=perEachOtherCase, **kwargs)
    
        an2.bootstrapDir = an2.workingDir + bootstrapDir

    if part1:
        if not df1main is None:
            an1.prepareDEG(df1main, df1other)

        an1.preparePerBatchCase(exprCutoff=exprCutoff1)
        an1.prepareBootstrapExperiments(parallel=parallelBootstrap)

        if perEachOtherCase:
            if not df2main is None:
                an2.prepareDEG(df2main, df2other)

            an2.preparePerBatchCase(exprCutoff=exprCutoff2)
            an2.prepareBootstrapExperiments(parallel=parallelBootstrap)

    if part2:
        an1.analyzeBootstrapExperiments()

        if perEachOtherCase:
            an2.analyzeBootstrapExperiments()
            an1.analyzeBootstrapExperiments()
            
        an1.reanalyzeMain()
        an1.analyzeCombinationVariant('Avg combo3avgs')
        an1.analyzeCombinationVariant('Avg combo4avgs')
        an1.bootstrapMaxpeakPlot('Avg combo3avgs')
        an1.bootstrapMaxpeakPlot('Avg combo4avgs')
        an1.analyzeAllPeaksOfCombinationVariant('Avg combo3avgs', nG=15, nE=15, fcutoff=0.5, width=50)
        an1.analyzeAllPeaksOfCombinationVariant('Avg combo4avgs', nG=15, nE=15, fcutoff=0.5, width=50)

        if doScramble:
            an1.scramble(['Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo3/', M=20)
            an1.scramble(['Markers', 'Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo4/', M=20)

        if perEachOtherCase:
            an2.reanalyzeMain()
            an2.analyzeCombinationVariant('Avg combo3avgs')
            an2.analyzeCombinationVariant('Avg combo4avgs')
            an2.bootstrapMaxpeakPlot('Avg combo3avgs')
            an2.bootstrapMaxpeakPlot('Avg combo4avgs')
            an2.analyzeAllPeaksOfCombinationVariant('Avg combo3avgs', nG=15, nE=15, fcutoff=0.5, width=50)
            an2.analyzeAllPeaksOfCombinationVariant('Avg combo4avgs', nG=15, nE=15, fcutoff=0.5, width=50)

            if doScramble:
                an2.scramble(['Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo3/', M=10)
                an2.scramble(['Markers', 'Binomial -log(pvalue)', 'Top50 overlap', 'Fraction'], subDir='combo4/', M=20)

    #an1.generateAnalysisReport()
    #
    #if perEachOtherCase:
    #    an2.generateAnalysisReport()

    if perEachOtherCase:
        return an1, an2

    return an1

if __name__ == '__main__':
  
    if platform.system()=="Windows":
        wdir = 'd:/Projects/A_Endothelial/VS/Endothelial/results/' 
    else:
        #wdir = '/mnt/research/piermarolab/Sergii/results/'
        wdir = '/mnt/scratch/PanglaoDCS1000/'
    
    if False:
        preparePerSampleDEgenes()
        prepareInput(wdir + 'PanglaoDB_byDCS_human/data.h5', 'Homo sapiens')
        prepareInput(wdir + 'PanglaoDB_byDCS_mouse/data.h5', 'Mus musculus')

    # Spearman: Human part1->1hr, mouse part1a->1hr&222GB, mouse part1b->10hrs&180GB
    parameters = dict(nCPUs=4 if platform.system()=="Windows" else 20, parallelBootstrap=False,
                PCNpath=os.path.join(os.path.dirname(__file__), 'data'), exprCutoff1=0.05, exprCutoff2=0.05, 
                genesOfInterest=receptorsListHugo_2555, knownRegulators=gEC22, perEachOtherCase=True, part1=False, part2=False,
                majorMetric='correlation',                # (1) correlation    (2) spearman 
                dendrogramMetric='euclidean',           # (1) euclidean    (2) correlation
                dendrogramLinkageMethod='ward')       # (1) ward    (2) complete (3) average

    parameters.update(dict(bootstrapDir='bootstrap%s/' % id))

    anHuman, anMouse = process_loc(*(None, None), *(None, None), wdir + 'PanglaoDB_byDCS_human_%s/' % parameters['majorMetric'], wdir + 'PanglaoDB_byDCS_mouse_%s/' % parameters['majorMetric'], **parameters)

    #anHuman.prepareBootstrapExperiments(parallel=False)
    anMouse.prepareBootstrapExperiments(parallel=False)
