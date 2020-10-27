.. _input-data:

**Input data format**
=====================

Expression data for **two different species** for comparison is required. For each of these species provide the input gene expression data is expected in one of the following formats:

1. Spreadsheet of comma-separated values ``csv`` where rows are genes, columns are cells with gene expression counts, this should be accompanied by another dataframe with two columns with one specifying batches and the other specifying corresponding cells.
Alternatively, the first row of the dataframe should be ``'batch'`` and the second ``'cell'``.

.. list-table:: 
        :header-rows: 0
        :align: left
        :widths: auto

        * - **Cell vs Genes**
          - **Batches and Cells**
        * - +-------+--------+--------+--------+--------+     
            | cell  | C1     | C2     | C3     | C4     |
            +=======+========+========+========+========+
            | G1    |        | 3      | 1      | 7      |
            +-------+--------+--------+--------+--------+
            | G2    | 2      | 2      |        | 2      |
            +-------+--------+--------+--------+--------+ 
            | G3    | 3      | 1      |        | 5      |
            +-------+--------+--------+--------+--------+
            | G4    | 10     |        | 5      | 4      |
            +-------+--------+--------+--------+--------+
            | ...   | ...    | ...    | ...    | ...    |
            +-------+--------+--------+--------+--------+
          - +-------+--------+
            | batch | cell   |
            +=======+========+
            | B1    | C1     | 
            +-------+--------+
            | B1    | C2     |
            +-------+--------+
            | B2    | C3     |
            +-------+--------+
            | B3    | C4     |
            +-------+--------+
            | ...   | ...    |
            +-------+--------+

or:

+-------+--------+--------+--------+--------+
| batch | batch0 | batch0 | batch1 | batch1 |
+-------+--------+--------+--------+--------+
| cell  | C1     | C2     | C3     | C4     |
+=======+========+========+========+========+
| G1    |        | 3      | 1      | 7      |
+-------+--------+--------+--------+--------+
| G2    | 2      | 2      |        | 2      |
+-------+--------+--------+--------+--------+
| G3    | 3      | 1      |        | 5      |
+-------+--------+--------+--------+--------+
| G4    | 10     |        | 5      | 4      |
+-------+--------+--------+--------+--------+
| ...   | ...    | ...    | ...    | ...    |
+-------+--------+--------+--------+--------+

2. ``Pandas DataFrame`` where ``axis 0`` is genes and ``axis 1`` are cells.
If the are batched in the data then the index of ``axis 1`` should have two levels, e.g. ``('batch', 'cell')``, 
with the first level indicating patient, batch or expreriment where that cell was sequenced, and the
second level containing cell barcodes for identification.

.. code:: python

    df = pd.DataFrame(data=[[2,np.nan],[3,8],[3,5],[np.nan,1]], 
                      index=['G1','G2','G3','G4'], 
                      columns=pd.MultiIndex.from_arrays([['batch0','batch1'],['C1','C2']], names=['batch', 'cell'])) 

