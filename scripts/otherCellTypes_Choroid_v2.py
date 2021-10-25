from decneo.commonFunctions import *
from decneo.analysisPipeline import process

if __name__ == '__main__':

    if False:
        try:
            # modes: ['ligands', 'receptors']
            # celltypes: ['Macrophage', 'Pericyte', 'SMC', 'Fibroblast', 'Endothelial']
            mode, celltype = sys.argv[1], sys.argv[2]
            print('mode/celltype:', mode, celltype, flush=True)
        except Exception as exception:
            print('ARG ERROR:', exception, flush=True)
            raise AttributeError

        workingDir = '/mnt/gs18/scratch/users/paterno1/otherCellTypes_choroid/'
        dataPath = processedDataDir + 'Voigt_all/' + 'df_Voigt_Choroid_AMD_and_Normal_expression_by_celltypes_{group}{celltype}.h5'
        PCN_path = '/mnt/ufs18/home-132/paterno1/Ben/'
        if mode=='ligands':
            gEC_choroidLigands =['ADIPOQ','ANGPT1','ANGPT2','ANGPTL3','BMP2','BMP7','C3','C4B','DLL1',
                     'FN1','IL13','IL4','JAG1','JAG2','LIF','OSM','S1P','SEMA3A','SEMA3C',
                     'SEMA3E','SEMA4A','SLIT2','TGFB1','TGFB2','TGFB3','VEGFA'] # 26
            knownRegulators = gEC_choroidLigands
            genesOfInterest = ligands_1777
            workingDirOther = 'otherCellTypes_lig/{}/Mus musculus/bootstrap/All/'
        elif mode=='receptors':
            knownRegulators = gEC22
            genesOfInterest = receptorsListHugo_2555
            workingDirOther = 'otherCellTypes_rec/{}/Mus musculus/bootstrap/All/'
        else:
            raise NotImplementedError

        an = process(pd.read_hdf(dataPath.format(group='', celltype=celltype), key='df'), 
                    pd.read_hdf(dataPath.format(group='non', celltype=celltype), key='df'), *(None, None),
                    workingDir + 'otherCellTypes_Choroid_AMD_and_Normal_LR_5_celltypes/%s/%s/' % (mode, celltype), 
                    workingDir + workingDirOther.format(celltype), 
                    nCPUs=20, parallelBootstrap=True, genesOfInterest=genesOfInterest, knownRegulators=knownRegulators, 
                    exprCutoff1=0.01, perEachOtherCase=False, nBootstrap=100, part1=True, part2=True, part3=True, PCNpath=PCN_path)

        #an.reanalyzeMain(togglePublicationFigure=False, includeClusterNumber=False, toggleIncludeHeatmap=False, toggleCalculateMeasures=False, toggleExportFigureData=True)
    if True:
        PCN_path = '/mnt/ufs18/home-132/paterno1/Ben/'
        dataPath = processedDataDir + 'Voigt_all/' + 'df_Voigt_Choroid_AMD_and_Normal_expression_by_celltypes_{group}{celltype}.h5'
        workingDir = '/mnt/gs18/scratch/users/paterno1/otherCellTypes_choroid/'
        for mode in ['ligands', 'receptors']:
            for celltype in ['Fibroblast', 'Endothelial','Macrophage', 'Pericyte', 'SMC']: #'Macrophage', 'Pericyte', 'SMC', 

                if platform.system() == "Windows":
                    workingDir = 'd:/Projects/A_Endothelial/VS/Endothelial/dev/otherCelltypes/'
                else:
                    workingDir = '/mnt/gs18/scratch/users/paterno1/otherCellTypes_choroid/'

                if mode=='ligands':
                    gEC_choroidLigands =['ADIPOQ','ANGPT1','ANGPT2','ANGPTL3','BMP2','BMP7','C3','C4B','DLL1',
                             'FN1','IL13','IL4','JAG1','JAG2','LIF','OSM','S1P','SEMA3A','SEMA3C',
                             'SEMA3E','SEMA4A','SLIT2','TGFB1','TGFB2','TGFB3','VEGFA'] # 26
                    knownRegulators = gEC_choroidLigands
                    genesOfInterest = ligands_1777
                    workingDirOther = 'otherCellTypes_lig/{}/Mus musculus/bootstrap/All/'
                elif mode=='receptors':
                    knownRegulators = gEC22
                    genesOfInterest = receptorsListHugo_2555
                    workingDirOther = 'otherCellTypes_rec/{}/Mus musculus/bootstrap/All/'
                else:

                    raise NotImplementedError

                an = process(pd.read_hdf(dataPath.format(group='', celltype=celltype), key='df'), 
                    pd.read_hdf(dataPath.format(group='non', celltype=celltype), key='df'), *(None, None),
                    workingDir + 'otherCellTypes_Choroid_AMD_and_Normal_LR_5_celltypes/%s/%s/' % (mode, celltype), 
                    workingDir + workingDirOther.format(celltype), 
                    nCPUs=20, parallelBootstrap=True, genesOfInterest=genesOfInterest, knownRegulators=knownRegulators, 
                    exprCutoff1=0.01, perEachOtherCase=False, nBootstrap=100, part1=True, part2=True, part3=True, PCNpath=PCN_path)
                an.reanalyzeMain(togglePublicationFigure=True, includeClusterNumber=False, toggleIncludeHeatmap=False, toggleCalculateMeasures=False, toggleExportFigureData=True)

    if False:
        PCN_path = '/mnt/ufs18/home-132/paterno1/Ben/'
        for mode in ['ligands', 'receptors']:
            for celltype in ['Fibroblast', 'Endothelial','Macrophage', 'Pericyte', 'SMC']: #'Macrophage', 'Pericyte', 'SMC', 

                if platform.system() == "Windows":
                    workingDir = 'd:/Projects/A_Endothelial/VS/Endothelial/dev/otherCelltypes/'
                else:
                    workingDir = '/mnt/gs18/scratch/users/paterno1/otherCellTypes_choroid/'

                if mode=='ligands':
                    gEC_choroidLigands =['ADIPOQ','ANGPT1','ANGPT2','ANGPTL3','BMP2','BMP7','C3','C4B','DLL1',
                             'FN1','IL13','IL4','JAG1','JAG2','LIF','OSM','S1P','SEMA3A','SEMA3C',
                             'SEMA3E','SEMA4A','SLIT2','TGFB1','TGFB2','TGFB3','VEGFA'] # 26
                    knownRegulators = gEC_choroidLigands
                    genesOfInterest = ligands_1777
                    workingDirOther = 'otherCellTypes_lig/{}/Mus musculus/bootstrap/All/'
                elif mode=='receptors':
                    knownRegulators = gEC22
                    genesOfInterest = receptorsListHugo_2555
                    workingDirOther = 'otherCellTypes_rec/{}/Mus musculus/bootstrap/All/'
                else:

                    raise NotImplementedError

                an = process(*(None, None), *(None, None),
                    workingDir + 'otherCellTypes_Choroid_AMD_and_Normal_LR_5_celltypes/%s/%s/' % (mode, celltype), 
                    workingDir + workingDirOther.format(celltype), 
                    nCPUs=1, parallelBootstrap=True, genesOfInterest=genesOfInterest, knownRegulators=knownRegulators, 
                    exprCutoff1=0.01, perEachOtherCase=False, nBootstrap=100, part1=True, part2=True, part3=True, PCNpath=PCN_path)
                an.reanalyzeMain(togglePublicationFigure=True, includeClusterNumber=False, toggleIncludeHeatmap=False, toggleCalculateMeasures=True, toggleExportFigureData=True)


                try:
                    an.analyzePerGeneCombinationVariant('Avg combo3avgs', hcutoff=0.2, fcutoff=0.1, width=50)
                    an.analyzePerGeneCombinationVariant('Avg combo3avgs', hcutoff=0.0, fcutoff=0.1, width=50)
                except:
                    pass
            
    if True:
        celltypes = ['Endothelial', 'Fibroblast', 'Pericyte', 'SMC', 'Macrophage']
        PCN_path = '/mnt/ufs18/home-132/paterno1/Ben/'
        for mode in ['receptors', 'ligands']:
            for celltype1 in celltypes:
                if platform.system() == "Windows":
                    workingDir = 'd:/Projects/A_Endothelial/VS/Endothelial/dev/otherCelltypes/'
                else:
                    workingDir = '/mnt/gs18/scratch/users/paterno1/otherCellTypes_choroid/'

                if mode=='ligands':
                    gEC_choroidLigands =['ADIPOQ','ANGPT1','ANGPT2','ANGPTL3','BMP2','BMP7','C3','C4B','DLL1',
                                'FN1','IL13','IL4','JAG1','JAG2','LIF','OSM','S1P','SEMA3A','SEMA3C',
                                'SEMA3E','SEMA4A','SLIT2','TGFB1','TGFB2','TGFB3','VEGFA'] # 26
                    knownRegulators = gEC_choroidLigands
                    genesOfInterest = ligands_1777
                    workingDirOther = 'otherCellTypes_lig/{}/Mus musculus/bootstrap/All/'
                elif mode=='receptors':
                    knownRegulators = gEC22
                    genesOfInterest = receptorsListHugo_2555
                    workingDirOther = 'otherCellTypes_rec/{}/Mus musculus/bootstrap/All/'
                else:
                    raise NotImplementedError

                panels = []
                for celltype2 in celltypes:
                    hcutoff, key = '0.0', 'choroid'
                    ddir = '/mnt/research/piermarolab/Sergii/LRpairs/v3_LR2422_Choroid100/'

                    if mode=='receptors':
                        a, b = celltype2, celltype1
                    else:
                        a, b = celltype1, celltype2

                    # a has ligands, b has receptors
                    tempName = ddir + 'v3_allpairs_%s_%s_%s_%s.h5' % (hcutoff, key, a, b)
                    data = pd.read_excel(tempName[:-3] + '_%s_max.xlsx' % mode[:-1], header=0, index_col=0)['count-1']

                    pname = '%s max' % celltype2
                    externalPanelsData.update({pname: data.loc[~data.index.duplicated(keep='first')]})
                    panels.append(pname)

                an = process(*(None, None), *(None, None),
                    workingDir + 'otherCellTypes_Choroid_AMD_and_Normal_LR_5_celltypes/%s/%s/' % (mode, celltype1), 
                    workingDir + workingDirOther.format(celltype1), 
                    nCPUs=1, parallelBootstrap=True, genesOfInterest=genesOfInterest, knownRegulators=knownRegulators, 
                    exprCutoff1=0.01, perEachOtherCase=False, nBootstrap=100, part1=False, part2=False, part3=False, PCNpath=PCN_path,
                    panels=['fraction', 'binomial', 'top50', 'markers', 'combo3avgs', 'combo4avgs'])

                an.reanalyzeMain(togglePublicationFigure=True, includeClusterNumber=False, toggleIncludeHeatmap=False, toggleCalculateMeasures=False, toggleExportFigureData=True)

    if True:
        cwd = '/mnt/gs18/scratch/users/paterno1/otherCellTypes_choroid/'
        #cwd = os.path.dirname(__file__).replace('\\', '/') + '/' # ''
        print(cwd)

        for celltype in ['Endothelial', 'Pericyte', 'Macrophage', 'SMC', 'Fibroblast']:
            for type in ['receptors','ligands']:
                for workingDir in ['otherCellTypes_Choroid_AMD_and_Normal_LR_5_celltypes/'+type + '/' + celltype + '/bootstrap/All/', 
                               'otherCellTypes_Choroid_AMD_and_Normal_LR_5_celltypes/'+ type + '/' + celltype + '/bootstrap/All/']:
                    for file in ['dendrogram-heatmap-correlation-data.xlsx', 'batches.txt']:
                    #for file in ['results All correlation euclidean ward.png', 'results All correlation euclidean ward.xlsx',
                    #            'all peaks Avg combo3avgs.png', 'All peaks Avg combo3avgs.xlsx', 'Avg combo3avgs_variant.xlsx',
                    #            'heatmap all peaks Avg combo3avgs.xlsx', 'near frequency Avg combo3avgs.xlsx']:
                    # 'dendrogram-heatmap-correlation-data.xlsx'
                        #source = cwd + workingDir + '%s/' % celltype
                        source = cwd + workingDir
                        destination = cwd + 'results 09 17 2021 0.0/' + type +'/'

                        if not os.path.exists(destination):
                            os.makedirs(destination)

                        if os.path.isfile(source + file):
                            shutil.copyfile(source + file, destination + '%s. ' % celltype + file)