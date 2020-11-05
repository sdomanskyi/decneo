# DECNEO

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

DECNEO is dependent on the following packages, so install or update: 
- [x] Matlibplot - for plotting results
- [x] NetworkX - for 
- [x] Pandas and NumPy and hdf5 - for processing data
- [x] rpy2 - for 
- [x] sklearn - for clustering 
- [x] tables - for clustering 

### Installation

To install DECNEO as a package:

	$ python setup.py install

Alternatively, you can create a local copy of this project for development purposes, and 
install the package from the cloned directory:

	git clone https://github.com/sdomanskyi/decneo
	python setup.py install

#### Troubleshooting 

NEEDS TO BE FILLED IN

### Deployment 

In your script import the package's main function ```process``` from the analysis pipeline:

	from decneo.analysisPipeline import process

Run the  ```process``` function. Here, for simplicity, we input our data and directories and use default values for the remaining parameters:

```python
process(pd.read_hdf(data, key='dfa'),
            pd.read_hdf(data, key='dfb'),
            None, None,                       
            'dir1/', 'dir2/')
```

During the initialization a number of parameters can be specified. For detailed list see documentation.

## Functionality 

### Overview

The main class, Analysis, includes tools for:

  1. **Pre-preprocessing**
  2. **Preparing bootstrap cases**
  3. **Calculating distance metric**
  4. **Generating dendrogram, heatmap, and panels for each case**
  4. **Analyzing combinations of measurements**
  5. **Identify peaks**
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

FILL IN

## Licensing 

DECNEO is released under an MIT License. Please also consult the folder LICENSES distributed with DECNEO regarding Licensing information for use of external associated content.

For development and use:
------------------------

> Note: this is not intended for public use yet. Publication is being prepared.

For development and testing of the documentation locally (on the development machine) install Sphinx by:

	$ pip install -U sphinx

To compile html version of the documentation:

	$ sphinx-build -E -a -b html ./docs/source ./docs/build

We are utilizing a 3rd party Sphinx extension sphinxcontrib-images extension, allowing to display documentation images in a organized way.