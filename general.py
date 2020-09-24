import os
import sys
import numpy as np
import pandas as pd
import h5py
import scipy.stats
import scipy.signal
import scipy.io
from scipy.spatial.distance import cdist, pdist, squareform
import time
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import plotly.graph_objects as go
import plotly.figure_factory as ff
import plotly.express as px
from plotly.offline import plot as plot_offline
from plotly.subplots import make_subplots
from sklearn.cluster import AgglomerativeClustering, KMeans, SpectralCoclustering
from sklearn.metrics import adjusted_rand_score, f1_score, roc_curve, auc, roc_auc_score
from scipy.cluster import hierarchy
from adjustText import adjust_text
from sknetwork.hierarchy import Paris, LouvainHierarchy
import copy
import networkx as nx
import multiprocessing
import platform
import warnings
import gzip
import pickle
import urllib.request
import tarfile
import json

from tables.exceptions import NaturalNameWarning
warnings.simplefilter("ignore", NaturalNameWarning)
#print('Ignoring any warnings of type: NaturalNameWarning')

if platform.system() == "Windows":
    RDataDirName = os.path.join('dev', 'PanglaoDBdata', '')
    MetadataDirName = os.path.join('dev', 'PanglaoDBdata', '')
    processedDataDir = os.path.join('data', 'pr', '')
else:
    RDataDirName = '/mnt/research/piermarolab/data/PanglaoDBh5/Raw/'
    MetadataDirName = os.path.join('')
    processedDataDir = '/mnt/research/piermarolab/Sergii/processedDataDir/'
    
def write(data, fileName, jsonFormat = False):
    
    '''Pickle object into a file.'''
    
    if jsonFormat:
        with gzip.GzipFile(fileName, 'w') as tempFile:
            tempFile.write(json.dumps(data).encode('utf-8'))

        return

    try:
        with gzip.open(fileName + '.pklz','wb') as tempFile:
            pickle.dump(data, tempFile, protocol=4)

    except Exception as exception:
        print(exception)

    return

def read(fileName, jsonFormat = False):

    '''Unpickle object from a file.'''

    if jsonFormat:
        with gzip.GzipFile(fileName, 'r') as tempFile:
            data = json.loads(tempFile.read().decode('utf-8'))

        return data

    if os.path.isfile(fileName + '.pklz'):
        try:
            with gzip.open(fileName + '.pklz','rb') as tempFile:
                data = pickle.load(tempFile)

                return data

        except Exception as exception:
            print(exception)

    return

