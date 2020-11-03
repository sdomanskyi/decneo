import pandas as pd
from decneo.analysisPipeline import process

demoData = 'demo/VoightChoroid4567RemappedData.h5'

if __name__ == '__main__':

    process(pd.read_hdf(demoData, key='dfa'),   # Endothelial cells
            pd.read_hdf(demoData, key='dfb'),   # Non-endothelial cells
            None, None,                         # Comparison dataset is provided
            'demo/',                            # Working directory
            'demo/fromPanglaoDBmouseAllbyDCS/', # Comparison dataset 
            parallelBootstrap=True,             # Set False if RAM is limited
            exprCutoff1=0.01,                   # Gene expression cutoff
            perEachOtherCase=False)             # Comparison mode setting