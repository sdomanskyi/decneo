import os
import gzip
import pickle
import urllib.request
import tarfile
import pandas as pd
import numpy as np
from scipy.io import mmread

from string import Template

def write(data, fileName):
    
    '''Pickle object into a file.'''

    try:
        with gzip.open(fileName + '.pklz','wb') as tempFile:
            pickle.dump(data, tempFile, protocol=4)

    except Exception as exception:
        print(exception)

    return

def read(fileName):

    '''Unpickle object from a file.'''

    if os.path.isfile(fileName + '.pklz'):
        try:
            with gzip.open(fileName + '.pklz','rb') as tempFile:
                data = pickle.load(tempFile)

                return data

        except Exception as exception:
            print(exception)

    return

def KeyInStore(key, file):

    try:
        with pd.HDFStore(file) as hdf5file:
            if "/" + key.strip("/") in hdf5file.keys():
                return True
            else:
                return False
    except Exception as exception:
        print(exception)

    return

def downloadFile(url, saveDir, saveName = None):

    if not os.path.exists(saveDir): 
        os.makedirs(saveDir)
    
    if saveName is None:
        saveName = url.strip('"').split('/')[-1:][0]

    path = os.path.join(saveDir, saveName)

    if os.path.isfile(path):
        print('File has been downloaded already')
    else:
        print('Downloading file:', url.strip('"'), end='\t', flush=True)

        try:
            urllib.request.urlretrieve(url.strip('"'), path)
            print('Done', flush=True)

        except Exception as exception:
            print(exception)

    return

def readRDataFile(fullPath, takeGeneSymbolOnly = True, saveToHDF = True, returnSizeOnly = False):

    if (not returnSizeOnly) and os.path.isfile(fullPath):
        print('File already exists:', fullPath)

        df = pd.read_hdf(fullPath, key='df')

        return df

    from rpy2.robjects import r as R

    R['load'](fullPath[:-len('.h5')])
    ls = np.array(R['ls']())
    #print('Variables in RData file:', ls)

    rSparseMatrix = R[ls[0]]

    print('Matrix size:', end='\t', flush=True)
    size = R('dim')(rSparseMatrix)
    print(size)

    if returnSizeOnly:

        return size[0], size[1]

    columns = pd.Index(np.array(R['colnames'](rSparseMatrix)))
    index = pd.Index(np.array(R['rownames'](rSparseMatrix)))

    if takeGeneSymbolOnly:
        index = index.str.split('_ENS', expand=True).get_level_values(0)
    
    R['writeMM'](obj=rSparseMatrix, file='%s.mtx' % (fullPath))

    df = pd.DataFrame(index=index, columns=columns, data=mmread('%s.mtx' % (fullPath)).toarray().astype(int))

    os.remove('%s.mtx' % (fullPath))

    if saveToHDF:
        df.to_hdf(fullPath, key='df', mode='a', complevel=4, complib='zlib')

        return df

    return df
