'''  Holds genes used in scRegulation and converts genes to the correct format
'''

from general import *

Mouse_to_Human_HUGO_conversion = read('geneLists/Mouse_to_Human_HUGO.gz', jsonFormat=True)

g4 = ['CD3D', 'CD28', 'CTLA4', 'PDCD1']

g16 = ['CD3D','CD27','CD28','CD82','ICOS','CD226','TNFRSF18','TNFRSF9', # stimulators
       'CTLA4','PDCD1','LAG3','BTLA','HAVCR2','VSIR','TIGIT','CD96']    # inhibitors

gEC23 = ['KDR','FLT1','FLT4','NRP1','NRP2','FGFR1','FGFR2','FGFR3','CXCR2','ROBO1',
         'ROBO4','ENG','PDGFRA','PDGFRB','TEK','KIT','MET','CLEC14A', # stimulators
         'CD36','CD47','VLDLR','PLXND1'] # inhibitors # on 08/24/2020 removed 'PNPLA2', due to low evidence that 'PNPLA2' is a receptor

g13 = ['PTPRB','TEK','TIE1','CD93','PECAM1','ENG','RAMP2','ADGRL4','CDH5','KDR','ESAM','ADGRF5','FLT1']

gAbove50_PanglaoMouse = ['ESAM','ADGRF5','RAMP2','FLT1','PECAM1','ENG','KDR','CD93','CDH5','ADGRL4','PTPRB','TIE1','TEK','CLEC14A','ROBO4']
gAbove50_PanglaoHuman = ['ESAM','RAMP2','PECAM1','CLEC14A','ADGRL4','ADGRF5','PTPRB','TIE1','FLT1','CDH5','ENG','CAV1','CALCRL', 
                        'CD93','KDR','TGFBR2','THBD','CD59','S1PR1','FCGRT','CD151','CD63','RAMP3']

pubMedAngiogenesisHits = {'KDR': 1853, 'FLT1': 681, 'PDGFRB': 546, 'NRP1': 237, 'TEK': 196, 'EPHB4': 173, 'ACVRL1': 153, 'FLT4': 120, 'EPHA2': 116, 'TIE1': 114, 'PECAM1': 100, 'NRP2': 100, 'APLNR': 100, 'ROBO4': 69, 'ACKR3': 69, 'S1PR1': 67, 'ROBO1': 67, 'LYVE1': 62, 'CDH5': 48, 'ADGRB1': 46, 'UNC5B': 43, 'FZD4': 41, 'CD248': 38, 'PLXND1': 31, 'PLXNB1': 29, 'PROKR1': 27, 'ITGAV': 23, 'CLEC14A': 22, 'PROKR2': 20, 'S1PR3': 20, 'CD93': 19, 'PLXNA1': 14, 'ADGRA2': 14, 'ADGRB3': 14, 'STAB1': 13, 'OGFR': 11, 'ADGRB2': 9, 'PTPRB': 7, 'S1PR5': 7, 'ADGRL4': 6, 'S1PR4': 4, 'GPR182': 3, 'PLXNB3': 3, 'PLXNC1': 3}

gECs = gEC23[:18]
gECi = gEC23[18:]

receptorsListHugo_2555 = np.loadtxt('geneLists/receptorsListHugo_2555.txt', dtype=str).tolist()

OSKM = ['POU5F1', 'SOX2', 'KLF4', 'MYC'] # OCT4, KLF4, SOX2, c-Myc
TFmarkers = OSKM + ['NANOG', 'GLIS1', 'NR5A2', 'SALL4'] # from 10.1152/physrev.00039.2017

if True:
    TF = np.loadtxt('geneLists/TF_1.01_HUGO.txt', dtype=str)
else:
    se = pd.read_excel('geneLists/TranscriptionFactors_DatabaseExtract_v_1.01.xlsx', index_col='HGNC symbol', header=0)['Is TF?']
    TF = np.unique(se[se == 'Yes'].index.values).tolist()
    import DigitalCellSorter
    Convert = DigitalCellSorter.DigitalCellSorter(matplotlibMode='TkAgg').gnc.Convert
    TF = np.unique(Convert(TF, 'alias', 'hugo', returnUnknownString=False))
    np.savetxt('geneLists/TF_1.01_HUGO.txt', TF, fmt='%s')


def populateExternalPanelsData(var):

    data = pd.read_excel('geneLists/Rate.xlsx', header=0, index_col=0)['Rate']
    var.update({'Evolutionary rate': data.loc[~data.index.duplicated(keep='first')]})

    data = pd.read_excel('geneLists/Age.xlsx', header=0, index_col=0)['Age']
    var.update({'Evolutionary age': data.loc[~data.index.duplicated(keep='first')]})

    data = pd.read_excel('geneLists/positive-negative-GO.xlsx', sheet_name='Lists', index_col=0, header=0, usecols=[0,1,2])
    var.update({'GOpositive': data['positive'].dropna().index.values.tolist()})
    var.update({'GOnegative': data['negative'].dropna().index.values.tolist()})

    var.update({'gAbove50_PanglaoMouse': gAbove50_PanglaoMouse})

    var.update({'gAbove50_PanglaoHuman': gAbove50_PanglaoHuman})

    var.update({'pubMed angiogenesis hits': pubMedAngiogenesisHits})

    return var

externalPanelsData = populateExternalPanelsData(dict())

choroid_paper_pan_markers = {
            'Schwann':          ['PLP1', 'MPZ', 'MBP', 'NCAM1', 'SCN7A', 'NGFR', 'L1CAM'],
            'Endothelial':      ['VWF', 'CD34', 'ICAM2'],
            'Melanocytes':      ['MLANA', 'PMEL', 'TYRP1', 'DCT'],
            'Pericytes/Smooth': ['RGS5'],
            'Fibroblasts':      ['FBLN1'],
            'RPE':              ['RPE65', 'BEST1'],
            'Lleukocytes':      ['PTPRC'],
            'B-cells':          ['CD79A'],
            'T-cells/NK':       ['CD2'],
            'Mast':             ['KIT'],
            'Macrophages':      ['AIF1']}

pan_markers_eye = {
    'Ganglion':             ['SNCG', 'POU4F1', 'NRN1', 'SLC17A6', 'ZIC1', 'RUNX1', 'MEF2C', 'IRX4', 'FST'],
    'Amacrine':             ['GAD1', 'SLC6A9', 'STX1B', 'CALB2', 'CHAT', 'TH', 'SLC17A8', 'EBF3', 'PROX1', 'PAX6', 'GAD1', 'SLC6A9'],
    'Schwann':              ['PLP1', 'MPZ', 'MBP', 'NCAM1', 'SCN7A', 'NGFR', 'L1CAM', 'SOX2', 'SOX10', 'POU3F1', 'GAP43', 'NGFR'],
    'Rod\n bipolar':        ['TRPM1', 'GRM6', 'SEBOX', 'PRKCA', 'VSTM2B', 'OTX2', 'VSX2'],
    'Cone\n bipolar':       ['SCGN', 'GRIK1', 'VSX1', 'LHX4', 'GLRA1', 'FEZF2', 'ZFHX4', 'OTX2', 'VSX2'],
    'Horiz.':               ['LHX1', 'CALB1', 'GJA10', 'ONECUT1', 'PROX1', 'PAX6'],
    'Cones':                ['OPN1SW', 'OPN1MW', 'GNAT2', 'ARR3', 'PDE6H'],
    'Rods':                 ['RHO', 'GNAT1', 'CNGA1', 'NRL', 'NR2E3'],
    'Astro.':               ['PAX2', 'GFAP', 'ALDH1L1', 'ALDOC', 'SLC1A2', 'S100B', 'AQP4'],
    'Muller\n glia':        ['SLC1A3', 'APOE', 'DKK3', 'GPR37', 'RLBP1', 'RAX', 'HES1', 'NOTCH1'],
    'RPE':                  ['RPE65', 'BEST1', 'RLBP1', 'TYR'],
    'Melano.':              ['MLANA', 'PMEL', 'TYRP1', 'DCT'],
    'Endo.':                ['PECAM1', 'CLDN5', 'CDH5', 'PECAM1', 'VWF', 'ERG', 'TEK', 'KDR', 'FLT1', 'CD34', 'ICAM2'],
    'Pericytes':            ['CSPG4', 'PDGFRB', 'ACTA1-', 'ACTA2-', 'DES-', 'MYL9-', 'MYH11-', 'GGT1', 'ABCC9', 'KCNJ8'],
    'Smooth\n muscle':      ['CSPG4', 'PDGFRB', 'ACTA1', 'ACTA2', 'DES', 'MYL9', 'MYH11', 'GGT1-', 'ABCC9-', 'KCNJ8-'],
    'Fibro.':               ['FBLN1', 'COL1A2', 'COL3A1', 'CXCL14', 'LRP1', 'PDPN', 'PROCR', 'S100A4', 'THY1', 'VEGFA', 'VEGFB', 'VTN'],
    'Leuk.':                ['PTPRC'],
    'B':                    ['CD79A', 'CD19', 'CD38'],
    'NK':                   ['CD2', 'CD3D-', 'CD3E-', 'CD3G-', 'FCGR3A', 'NCAM1'],
    'T':                    ['CD2', 'CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B'],
    'Macr.':                ['AIF1', 'ITGAM', 'CD68', 'CD163'],
    'mDC':                  ['CD1C', 'CD83', 'THBD', 'CD209'],
    }

Heng_eye_pan_markers = {
    'Amacrine':             ['GAD1', 'SLC6A9', 'STX1B', 'CALB2', 'CHAT', 'TH', 'SLC17A8', 'EBF3'],
    'Astrocytes':           ['PAX2', 'GFAP', 'GLUL'],
    'Cone bipolar':         ['SCGN', 'GRIK1', 'VSX1', 'LHX4', 'GLRA1', 'FEZF2', 'ZFHX4', 'OTX2', 'VSX2'],
    'Cones':                ['OPN1SW', 'OPN1MW', 'GNAT2', 'ARR3', 'PDE6H'],
    'Horizontal':           ['LHX1', 'CALB1', 'GJA10', 'ONECUT1'],
    'Muller glia':          ['SLC1A3', 'APOE', 'DKK3', 'GPR37', 'RLBP1', 'RAX', 'HES1', 'NOTCH1', 'GLUL'],
    'Perivascular':         ['MYL9', 'CSPG4', 'PDGFRB', 'MYH11', 'DES', 'ACTA2'],
    'Retinal ganglion':     ['SNCG', 'POU4F1', 'NRN1', 'SLC17A6'],
    'Rod bipolar':          ['TRPM1', 'GRM6', 'SEBOX', 'PRKCA', 'VSTM2B'],
    'Rods':                 ['RHO', 'GNAT1', 'CNGA1', 'NRL', 'NR2E3'],
    'Vascular endothelial': ['CLDN5', 'CDH5', 'PECAM1', 'VWF', 'ERG', 'TEK', 'KDR', 'FLT1']}


pan_names_eye_dict = {
    'Astrocytes': 'Astro.', 
     'Microglia': 'Astro.', 
     'Muller glia': 'Muller\n glia', 
     'Schwann cells': 'Schwann', 
     'Amacrine cells': 'Amacrine', 
     'Ganglion cells': 'Ganglion', 
     'Cone bipolar cells': 'Cone\n bipolar', 
     'Rod bipolar cells': 'Rod\n bipolar', 
     'Horizontal cells': 'Horiz.', 
     'Cones': 'Cones', 
     'Rods': 'Rods', 
     'RPE cells': 'RPE', 
     'Melanocytes': 'Melano.', 
     'Fibroblasts': 'Fibro.', 
     'Vascular endothelial cells': 'Endo.', 
     'Pericytes': 'Smooth\n muscle', 
     'Smooth muscle cells': 'Smooth\n muscle', 
     'T cells': 'T', 
     'B cells': 'B', 
     'NK cells': 'NK', 
     'mDCs': 'Macr.', 
     'pDCs': 'Macr.', 
     'Monocytes': 'Macr.', 
     'Macrophages': 'Macr.', 
 }

pan_names_eye_list = [
 'Ganglion',
 'Amacrine',
 'Schwann',
 'Rod\n bipolar',
 'Cone\n bipolar',
 'Horiz.',
 'Cones',
 'Rods',
 'Astro.',
 'Muller\n glia',
 'RPE',
 'Melano.',
 'Endo.',
 'Pericytes',
 'Smooth\n muscle',
 'Fibro.',
 'Leuk.',
 'B',
 'NK',
 'T',
 'Macr.',
 'mDC',
 ]

