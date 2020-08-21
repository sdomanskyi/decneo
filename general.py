import os
import sys
import numpy as np
import pandas as pd
from io import *
from PanglaoDBannotation import *
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
import DigitalCellSorter
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

from tables.exceptions import NaturalNameWarning
warnings.simplefilter("ignore", NaturalNameWarning)
print('Ignoring any warnings of type: NaturalNameWarning')

if platform.system() == "Windows":
    RDataDirName = os.path.join('dev', 'PanglaoDBdata', '')
    MetadataDirName = os.path.join('dev', 'PanglaoDBdata', '')
    processedDataDir = os.path.join('')
else:
    RDataDirName = '/mnt/research/piermarolab/data/PanglaoDBh5/Raw/'
    MetadataDirName = os.path.join('')
    processedDataDir = '/mnt/research/piermarolab/Sergii/processedDataDir/'
