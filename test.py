from commonFunctions import *
from analysisPipeline import Analysis

if __name__ == '__main__':

    print('pandas:', pd.__version__)

    an = Analysis(genesOfInterest=receptorsListHugo_2555, knownRegulators=gEC23, panels=['combo3avgs', 'combo4avgs', 'fraction', 'binomial', 'markers', 'top50'], workingDir='results/PanglaoDB_lung_mouse/', otherCaseDir='/mnt/research/piermarolab/Sergii/Endothelial by PanglaoDB definition/EC all bootstrap 100 w21/Homo sapiens/All/')

    an.compareTwoCases(an.bootstrapDir + 'All/', an.otherCaseDir, name1='name1', name2='name2', saveName=os.path.join(an.bootstrapDir + 'All/comparison'))
    an.reanalyzeMain()
