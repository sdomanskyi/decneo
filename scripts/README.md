# Supplementary and Usage Example Scripts

This folder contains scripts that show how to use DECNEO. 

## Supplement 

``PanglaoDBannotateWithDCS.py`` utilizes `Digital Cell Sorter` to annotate the Panglao DB cells for input for DECNEO. 

> DCS can be found here: https://github.com/sdomanskyi/DigitalCellSorter

``SRAs_Alona_DCS_analysis.py`` TEXT

## Usage  

### Demo 

We have made an example execution file ```demo.py``` that shows how to use DECNEO.

Download file ``VoightChoroid4567RemappedData.h5`` (456.7 Mb) 
from https://figshare.com/DECNEO/

> Processed by ``processRemappedChoroid.py``, this file contains normalized gene expression of 
> 27504 genes of 7996 endothelial cells from 8 batches, and 5704 non-endothelial cells from 8 batches. 
> Genes that are not expressed in endothelial cells are removed from non-endothelial cells dataset

Save the downloaded data file to ``demo/``, or otherwise modify path in ``demoData`` of
``demo.py``:

To execute the complete script ```demo.py``` run:

	python demo.py
  
Take a closer look at ``demo.py`` to see how DECNEO can be run.  

### Other Examples 

The folder also contains many other example scripts where we run DECNEO on different cell types of the Panglao DB. 
There are also available examples showing how to analyze only specified SRAs. 


