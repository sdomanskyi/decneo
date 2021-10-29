# Supplementary and Usage Example Scripts

This directory contains scripts that use DECNEO. 

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

The directory also contains many other scripts where we run DECNEO on different cell types of the PanglaoDB. 
There are also available scripts showing how to analyze only specified SRAs. 

### Ramilowski Pair Heatmap
These scripts generate the 2d heatmaps of ligand-receptor Ramilowski pairs.

1. Dendrograms of each cell type are created. This is done first on the Panglao samples with ``otherCellTypes_v2.py`` and then on the Choroid Samples with ``otherCellTypes_Choroid_v2.py``.
2. Ramilowski (RL) pairs of ligand-receptors pairs for each cell type can be created. The output is an excel file containing a specific pair and its neighboring pairs in the dendrogram.
3. ``make_sheet_from_all_pairs.py`` finds the number of unique RL pairs in each of the neighboring pairs. A unique RL pair is defined as a neighboring RL pair that does not contain the ligand or receptor of the specific pair.
4. ``heatmap.py`` contains functions that create heatmap of a ligand-receptor RL pair. A heatmap of ligand dendrogram values against receptor dendrogram values can be seen with the function ``heatmap``.
Another interesting function is ``heatmap_avg3sum`` which takes the average of the ligand dendrogram values and receptor dendrogram values (in the figure from the paper this is 'Avg Combination of measures'). 
This value is scaled from 0 to 1. The correlation of this average against the number of unique RL pairs can be visualized with the function ``Unique_avg3sum_Corr``. There is also a mode to create a heatmap from products of RL pairs.
``MakeBounds`` function group neighboring non-zero values in the heatmap. On default, the box containing the maximum value in the heatmap is only displayed.
``Avg3andUnique_Heatmap`` creates a heatmap via the average of RL pairs and scaled unique RL pairs. Scaled unique RL pairs is the number of unique RL pairs/maximum number of unique RL pairs in a specific RL cell type. Since both quantities are scaled from 0 to 1, then the average of these two values also range from 0 to 1. ``MakeNoNa_RLavg3`` makes any average less or equal to a specified threshold (0.5 by default) 0.
5. ``DBSCAN_fig.py`` creates a similar heatmap with ``DBSCAN`` from ``sklearn``. Non-zero values are colored and grouped based on specific clusters. On default, the cluster containing the largest value in the heatmap is only shown.
6. ``newfig.py`` creates the combined figure seen in the paper. By default, the ligand dendrogram is on the y-axis and the receptor dendrogram is on the x-axis. Also, the bounds from the DBSCAN figure are drawn on this heatmap. This script also does the necessary operations from ``heatmap.py`` and ``DBSCAN_fig.py``.

To generate a similar figure shown in the paper, the following should be run in order:
1. ``otherCellTypes_v2.py``
2. ``otherCellTypes_Choroid_v2.py``
3. ``make_sheet_from_all_pairs.py``
4. ``newfig.py``





