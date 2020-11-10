# DECNEO

> :warning: Warning: this software is not intended for public use yet. Publication is being prepared.

This repository contains DECNEO, a Python package that provides bioinformatics utilities for analyzing single cell transcriptomics datasets. DECNEO implements in silico detection of transcriptional regulation genes. The documentation is available at Read the Docs: https://decneo.readthedocs.io/

![logo](https://github.com/sdomanskyi/decneo/blob/master/docs/source/DECNEO.svg)

- [Getting Started](#getting-started)
  * [Dependencies](#dependencies)
  * [Installation](#installation)
    + [Troubleshooting](#troubleshooting)
  * [Deployment](#deployment)
- [Functionality](#functionality)
  * [Overview](#overview)
  * [Input Data Format](#input-data-format)
  * [Usage Example](#usage-examply)
  * [Output](#output)
- [Funding](#funding)
- [Licensing](#licensing)

## Getting Started

These instructions will get you a copy of the project up and running on your machine for data analysis, development or testing purposes.

### Dependencies 

The software runs in Python >= 3.8

DECNEO is dependent on the following packages, that are installed/updated with installation of DECNEO: 
- [x] Matplotlib - plotting from Python
- [x] NetworkX - used in network enrichment analysis
- [x] Pandas and tables - for data storage and analysis
- [x] NumPy - for processing data
- [x] sklearn - we use clustering algorithms and metrics
- [x] adjustText - optimization of text labels locations in plots

### Installation

To install DECNEO as a package:

	$ pip install decneo

Alternatively, you can clone a local copy of this project, and 
install the package from the cloned directory:

	git clone https://github.com/sdomanskyi/decneo
	python setup.py install

#### Troubleshooting 

If reading demo data gives error "unsupported pickle protocol: 5" make sure that python 3.8 is used and 
latest version of pandas and tables is installed.

## Functionality 

### Overview

The main class, Analysis, includes tools for:

  1. **Pre-preprocessing**
  2. **Preparing bootstrap cases**
  3. **Calculating distance metric**
  4. **Generating dendrogram, heatmap, and panels for each case**
  4. **Analyzing combinations of measurements**
  5. **Identifying peaks**
  6. **Recording peak genes**

### Input Data Format 

Expression data for **two different** species for comparison is required. 
For each of these species provide the input gene expression data is expected in one of the following formats:

1. Spreadsheet of comma-separated values ``csv`` where rows are genes, columns are cells with gene expression counts, this should be accompanied by another dataframe with two columns with one specifying batches and the other specifying corresponding cells.
Alternatively, the first row of the dataframe should be ``'batch'`` and the second ``'cell'``.

2. ``Pandas DataFrame`` where ``axis 0`` is genes and ``axis 1`` are cells.
If the are batched in the data then the index of ``axis 1`` should have two levels, e.g. ``('batch', 'cell')``, 
with the first level indicating patient, batch or expreriment where that cell was sequenced, and the
second level containing cell barcodes for identification.

For examples refer to documentation. 

### Usage Example 

We have made an example execution file ```demo.py``` that shows how to use ```decneo``.

Download file ``VoightChoroid4567RemappedData.h5`` (456.7 Mb) 
from https://figshare.com/DECNEO/

> This file contains normalized gene expression of 27504 genes of 7996 endothelial cells from 
> 8 batches, and 5704 non-endothelial cells from 8 batches. Genes that are not expressed in 
> endothelial cells are removed from non-endothelial cells dataset

Save the downloaded data file to ``demo/``, or otherwise modify path in ``demoData`` of
``demo.py``:

See details of the script ```demo.py``` at:

> [Example walkthrough of demo.py script](https://github.com/sdomanskyi/decneo/blob/master/scripts/demo.py)

To execute the complete script ```demo.py``` run:

	python demo.py

### Output 

Outputs all resulting directories, files, and figures to directory specified as the ``workingDir`` when creating an instance of class ``Analysis``. 
It will also output an analysis report detailing all results and figures.

For a detailed list, refer to the documentation. 

## Funding 

This research project is a part of R01GM122085 grant, funded by NIH/NIGMS.

## Licensing 

DECNEO is released under an MIT License. Please also consult the folder LICENSES distributed with DECNEO regarding Licensing information for use of external associated content.
