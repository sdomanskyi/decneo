''' Common functions used in analysis pipeline 
'''

from .general import *
from .genes import *

cleanListString = lambda c: str(list(c)).replace(' ', '').replace("'", '').replace(']', '').replace('[', '').replace(',', ', ')

def movingAverageCentered(a, halfWindowSize, looped = False):

    '''Function to smooth a 1d signal
        
    Parameters:
        a: ndarray
            Input data

        halfWindowSize: int
            Size of half-window for averaging

        looped: boolean, Default False
            Determined looped behaviour at the boundaries

    Returns:
        ndarray
            Smoothed signal

    Usage:
        movingAverageCentered(a, halfWindowSize)
    '''

    '''
    An intuitive way:
    b = np.append(np.append(a[-n:], a), a[:n])
    a_smoothed = np.array([np.sum(b[i-n:i+n+1]) / (2.*n + 1.) for i in range(n, len(b)-n)])

    Faster way is to use np.cumsum(), especially for long vectors

    if False:
        np.random.seed(0)
        a = np.sin(np.random.rand(100)*2.*180./np.pi) + 1.
        a[a > 0.5] = 0
        print(np.sum(a))
        plt.plot(np.arange(len(a)), a, 'o')

        for n in [1,2,5,10]:
            a_avg = movingAverageCenteredLooped(a, n, looped=False)
            print(np.sum(a_avg))
            plt.plot(np.arange(len(a)), a_avg)

        plt.show()
    '''

    n = halfWindowSize

    if looped:
        b = np.append(np.append(a[-n:], a), a[:n])
    else:
        b = np.append(np.append(np.zeros(n), a), np.zeros(n))

    s = np.cumsum(b) / (2.*n + 1.)
    s[2*n+1:] -= s[:-2*n-1]

    aSm = s[2*n:]

    if not looped:
        for i in range(n):
            e = i + n + 1
            aSm[i] = a[:e].mean()
            aSm[len(a) - i - 1] = a[-e:].mean()

    return aSm

def get_mean_std_cov_ofDataframe(df):

    '''Calculate mean, standard deviation, and covariance 
        
    Parameters:
        df: pandas.DataFrame
            Data with bootstrap experiment data 
        
    Returns:
        pandas.DataFrame:
            DataFrame with columns expressing mean, standard deviation, and covariance of input columns

    Usage:
        get_mean_std_cov_ofDataFrame(df)
    '''

    se_mean = df.mean(axis=1)
    se_std = df.std(axis=1)

    df['1 mean'] = se_mean
    df['2 std'] = se_std
    df['3 cov'] = se_std / se_mean
    df = df[['1 mean', '2 std', '3 cov']]

    df.columns = ['mean', 'std', 'cov']

    df = df.sort_values('mean', axis=0, ascending=False)

    return df

def getGenesOfPeak(se, peak=None, heightCutoff = 0.5, maxDistance = None):

    '''Find peak region of greatest value
        
    Parameters:
        se: Series
            Normalized aggregated data 

        peak: ndarray, Default None
            Indices of max values in data

        heighCutoff: float, Default 0.5
            Height/value considered to be in peak 

        maxDistance: int, Default None
            Maximum distance away considered to be in peak
        
    Returns:
        Genes appearing the peak 

    Usage:
        getGenesOfPeak(se)
    '''

    se /= se.max()

    if peak is None:
        i0 = np.argmax(se.values)
    else:
        i0 = peak

    genes = [se.index[i0]]

    if heightCutoff is None:
        heightCutoff = -np.inf

    if maxDistance is None:
        r = range(i0 + 1, i0 + 2*10**3, 1)
    else:
        r = range(i0 + 1, i0 + 1 + maxDistance, 1)
    for i in r:
        try:
            if se.iloc[i] >= (heightCutoff * se.iloc[peak]):
                genes.append(se.index[i])
            else:
                break
        except:
            break

    if maxDistance is None:
        r = range(i0 - 1, 0, -1)
    else:
        r = range(i0 - 1, max(0, i0 - 1 - maxDistance), -1)
    for i in r:
        try:
            if se.iloc[i] >= (heightCutoff * se.iloc[peak]):
                genes.append(se.index[i])
            else:
                break
        except:
            break

    return genes

def getPeaks(se, threshold = 0.2, distance = 50, prominence = 0.05):

    '''Find peak regions 
        
    Parameters:
        se: Series
            Normalized aggregated data 

        threshold: float, Default 0.2
            Minimum value to be considered peak 

        distance: int, Default 50
            Minimum horizontal distance between peaks 

    Returns:
        Indices of peaks satisfying input conditions 

    Usage:
        getPeaks(se)
    '''

    se /= se.max()

    peaks = scipy.signal.find_peaks(np.hstack([0, se.values, 0]), distance=distance, prominence=prominence)[0] - 1
    peaks = peaks[se[peaks] >= threshold]

    return peaks

def getDistanceOfBatch(args):

    '''Calculate correlation distance of metric for given batch 

    Parameters:
        args: tuple
            Tuple that contains: 
                batch: str
                    Batch identifier
                    
                df_sample: pandas.DataFrame
                    Expression data 

                metric: str
                    Metric name (e.g. 'correlation') 

                genes: list or 1d numpy.array
                    Genes of interest

                minSize: int
                    Minimum size of input pandas.DataFrame

                cutoff: float
                    Cutoff for percent expression of input data 

    Returns:
        tuple:
            Results in form of a tuple:
                pandas.Series 
                    Series containting correlation distance 

                str
                    Batch identifier 

                pandas.Series
                    Series of genes 

    Usage:
        getDistanceOfBatch
    '''

    batch, df_sample, metric, genes, minSize, cutoff = args

    try:
        df_sample = df_sample.loc[((df_sample > 0).sum(axis=1) / df_sample.shape[1]) >= cutoff]
    
        if df_sample.shape[1] < minSize:

            return pd.DataFrame(), '', set()

        print('\t', batch, df_sample.shape, end=' ')

        temp_genes = df_sample.index.intersection(genes)
        print(len(temp_genes), end=' ', flush=True)

        measure = cdist(df_sample.loc[temp_genes].values, df_sample.values, metric=metric).T

        se_batch_corr = pd.DataFrame(measure, index=df_sample.index, columns=temp_genes).stack().reorder_levels([1, 0]).sort_index()
        se_batch_corr = pd.concat([se_batch_corr], keys=[batch], axis=1, sort=False)
        se_batch_corr = pd.concat([se_batch_corr], keys=[metric], axis=1, sort=False)

    except Exception as exception:
        print('\nERROR in correlation calculation:', exception, '\n')

        return pd.DataFrame(), '', set()

    return se_batch_corr.copy(), batch, pd.Series(index=se_batch_corr.index)

def reindexSeries(args):

    '''Assists in reindexing Series 

    Parameters:
        args: tuple
            Tuple that contains: 
                se: pandas.Series
                    Series to perform reindexing 
                    
                batch: str
                    Batch identifier 

                index: list or 1d numpy.array
                    List of genes 

    Returns:
        tuple:
            Results in form of a tuple:
                pandas.Series
                    Reindexed pandas.Series

                str
                    Batch identifier


    Usage:
        reindexSeries
    '''

    se, batch, index = args

    print('\t', batch, end=' ', flush=True)
    se = se.replace(0., np.nan).reindex(index)

    return se, batch

def get_df_distance(df, metric = 'correlation', genes = [], analyzeBy = 'batch', minSize = 10, groupBatches = True, pname = None, cutoff = 0.05, nCPUs = 4):

    '''Calculate distance measurement
        
    Parameters:
        df: pandas.DataFrame
            Expression data 

        metric: str, Default 'correlation'
            Metric name (e.g. 'correlation') 

        genes: list, Default []
            Genes for analysis

        analyzeBy: str, Default 'batch'
            Level to analyze data by (e.g. batches)

        minSize: int, Default 10
            Minimum size of input pandas.DataFrame

        groupBatches: boolean, Default True 
            Whether to group batched or save per-batch distance measure 

        pname: Default None
            Deprecated

        cutoff: float, Default 0.05
            Cutoff for percent expression of input data 

        nCPUs: int, Default 4
            Number of CPUs to use for multiprocessing

    Returns:
        pandas.DataFrame 
            Distance measure

    Usage:
        get_df_distance(df)
    '''

    print('\tMetric:', metric, '\tAnalyzeBy:', analyzeBy, '\tminSize:', minSize)
    print('df_expr', '%1.1fMb'%(sys.getsizeof(df) / 1024**2), flush=True)

    genes = np.unique(genes)
    print('\tReceived %s genes for analysis' % len(genes))

    temp_batches = np.unique(df.columns.get_level_values(analyzeBy).values)
    print('Staring multiprocessing for %s batches' % len(temp_batches), flush=True)
    pool = multiprocessing.Pool(processes=nCPUs)
    results = pool.map(getDistanceOfBatch, [(batch, df.xs(key=batch, level=analyzeBy, axis=1, drop_level=False), metric, genes, minSize, cutoff) for batch in temp_batches])
    pool.close()
    pool.join()
    print()

    del df

    index = pd.concat([item[2] for item in results if item[1] != ''], axis=0).index.drop_duplicates()
    print('\tIndex pairs:', len(index), '\tSelected genes:', len(index.levels[0]), '\tAll genes:', len(index.levels[1]))

    print('Starting multiprocessing for reindexing', flush=True)
    pool = multiprocessing.Pool(processes=nCPUs)
    results = pool.map(reindexSeries, [(item[0], item[1], index) for item in results if item[1] != ''])
    pool.close()
    pool.join()
    print()

    print('Merging reindexed batches', flush=True)
    df_batches = pd.DataFrame(index=index, dtype=float)
    for se, batch in results:
        print('\t', batch, end=' ', flush=True)
        df_batches[batch] = se
    print()

    del results
    
    print(df_batches, flush=True)
    print('df_batches', '%1.1fMb'%(sys.getsizeof(df_batches) / 1024**2), flush=True)

    if groupBatches:
        df_batches = pd.Series(data=np.nanmedian(df_batches.values, axis=1, overwrite_input=True), index=df_batches.index).unstack(0)

        for gene in df_batches.columns:
            try:
                df_batches.loc[gene, gene] = 0.
            except:
                pass

    return df_batches

def reduce(vIn, size = 100):

    '''Interpolate data to reduce the number of data points
        
    Parameters:
        vIn: 1d vector
            Data to resample

        size: int, Default 100
            New data size
        
    Returns:
        Resampled data

    Usage:
        reduce(vIn)
    '''

    '''
    def reduce(v, size = 100):

        bins =  np.linspace(np.min(v), np.max(v), num=size)

        return bins[np.digitize(v, bins) - 1]
    '''

    v = vIn.copy()

    wh = np.where(~np.isnan(v))[0]

    if len(wh) > 0:
        bins =  np.linspace(np.min(v[wh]), np.max(v[wh]), num=size)

        v[wh] = bins[np.digitize(v[wh], bins) - 1]
    else:
        v[0] = 0.

    return v

def metric_euclidean_missing(u, v):

    '''Metric of euclidean distance between two arrays, excluding missing points
        
    Parameters:
        u: 1d vector
            Data array

        v: 1d vector
            Data array
        
    Returns:
        ndarray 
            Non-negative squareroot of the array, element-wise

    Usage:
        metric_euclidean_missing(u, v)
    '''

    wh = np.where(~np.isnan(u * v))[0]
                    
    return np.sqrt(((u[wh] - v[wh])**2).sum())

def binomialEnrichmentProbability(nx_obj, enriched_genes, target_genes = False, background_genes = False, PCNpath = 'data/'):

    '''Takes in a network x object and list of enriched genes and calculates the binomial
    enrichment based on the number of enriched interaction there are.
    Uses Survival function (also defined as 1 - cdf): scipy.stats.binom.sf(k, n, p, loc=0)

    Parameters:
        nx_obj: networkx.Graph or str
            A networkx undirected graph or edge list file name.

        enriched_genes: list
            List of enriched genes.

        target_genes: list or boolean, Default False
            List of target_genes. Default use all genes in background.

        background_genes: list or boolean, Default False
            List of genes to use as background probability. Default use all genes in the network.

    Returns:
        pandas.DataFrame
    '''
    
    if not isinstance(nx_obj, nx.Graph):
        if isinstance(nx_obj, str):
            if not os.path.isfile(os.path.join(PCNpath, 'PCN.pklz')):
                try:
                    nx_obj = nx.read_edgelist(nx_obj).to_undirected()
                except Exception as exception:
                    print(exception)

                write(nx_obj, os.path.join(PCNpath, 'PCN'))
            else:
                nx_obj = read(os.path.join(PCNpath, 'PCN'))
        else:
            print("Error: nx_obj needs to be a networkx graph or a network edgelist file name")
            return

    if background_genes is False:
        background_genes = nx_obj.nodes()
        
    if target_genes is False:
        target_genes = background_genes

    enriched_genes = set(enriched_genes).intersection(background_genes)
    
    enr_with_enr = 0
    enr_non_enr = 0
    for gene in enriched_genes:
        gene_neighbors = set(nx_obj.neighbors(gene))
        enr_non_enr += len(gene_neighbors - enriched_genes)
        enr_with_enr += len(gene_neighbors.intersection(enriched_genes))
        
    # Divide by 2 to not double count enriched interaction
    enr_with_enr = int(enr_with_enr/2)
    
    p = float(enr_non_enr + enr_with_enr)/float(len(nx_obj.edges())) 
    nx_nodes = set(nx_obj.nodes())
    df_binom = pd.DataFrame()
    for gene in target_genes:
    
        if gene not in nx_nodes: 
            continue

        gene_neighbors = set(nx_obj.neighbors(gene))
        num_neighbors_enriched = len(gene_neighbors.intersection(enriched_genes))
        num_neighbors_background = len(gene_neighbors.intersection(background_genes))

        if num_neighbors_background == 0: 
            continue
        
        df_binom.loc[gene, 'Enriched_Interact'] = num_neighbors_enriched
        df_binom.loc[gene, 'Total_Interact'] = num_neighbors_background
        df_binom.loc[gene, 'Binomial_Prob'] = scipy.stats.binom.sf(num_neighbors_enriched, num_neighbors_background, p, loc=1)
        
    return df_binom

def getROC(data):

    '''Calculate axes of ROC (false positive rates and true positive rates) using index as thresholds
        
    Parameters:
        data: 1d vector
            Input data
        
    Returns:
        np.array 
            Array holding false positive rates

        np.array
            Array holding true positive rates

    Usage:
        getROC(data)
    '''

    fpr = []
    for i in range(len(data)):
        FP = np.sum(data[:i] == False)
        TN = np.sum(data[i:] == False)
        fpr.append(FP/(FP+TN))

    tpr = []
    for i in range(len(data)):
        TP = np.sum(data[:i] == True)
        FN = np.sum(data[i:] == True)
        tpr.append(TP/(TP+FN))

    return np.array(fpr), np.array(tpr)

def clusterData(data, n_clusters = None, method = 'Spectral', random_state = None):

    '''Return cluster labels
    
    Parameters: 
        data: np.array 
            Array of shape (features, objects)

        n_clusters: int, Default None 
            Number of desired clusters 

        method: str or int , Default 'Spectral'
            Method used cluster the data 
            Options: 1 or 'Agglomerative', 2 or 'Spectral', 3 or 'KMeans'

        random_state: int, Default None
            Used to determine randomness deterministic 
            Methods 2 and 3 initial state is random, unless random_state is specified

    Returns:
        list
            Cluster assignment of each object

    Usage:
        clusterData(data)

    '''

    if n_clusters is None:
        
        print('Specify number of clusters via n_clusters')
        return

    if method in ['Agglomerative', 1]:
        model = AgglomerativeClustering(n_clusters=n_clusters).fit(data.T)
        labels = model.labels_

    elif method in ['Spectral', 2]:
        model = SpectralCoclustering(n_clusters=n_clusters, random_state=random_state)
        model.fit(data.T)
        labels =  model.row_labels_

    elif method in ['KMeans', 3]:
        model = KMeans(n_clusters=n_clusters, random_state=random_state).fit(data.T)
        labels =  model.labels_

    print('Cluster sizes:', np.unique(labels, return_counts=True)[1])
                            
    return labels

def testNumbersOfClusters(data, text = '', n_min = 2, n_max = 20, k = 10):

    '''Test different number of clusters with KMeans, calculate ARI for each
    
    Parameters: 
        data: np.array 
            Array of shape (features, objects)

        text: str, Default ''
            String identifier

        n_min: int, Default 2
            Minimum number of cluster

        n_max: int, Default 20 
            Maximum number of cluster

        k: int, Default 10 
            Number of iterations for each clusters number

    Usage: 
        testNumbersofClusters(data)
    '''

    print()

    for n in range(n_min, n_max):
        labels = [KMeans(n_clusters=n).fit(data.T).labels_ for i in range(k)]

        score = np.zeros((k, k))
        for i in range(k):
            for j in range(k):
                score[i, j] = adjusted_rand_score(labels[i], labels[j])

        score[np.triu_indices(k, k=0, m=k)] = np.nan
        score = np.nanmean(score.flatten())

        print('%s, n=%s, k=%s, ARI=%s' % (text, n, k, score), flush=True)

    print()

    return

def KeyInStore(key, file):

    '''Check if a key present in HDF store file

    Parameters:
        key: str
            The key to check

        file: str
            Path to HDF store file

    Returns
        True or False
            Whether key is in store or not
    
    Usage:
        KeyInStore(key, file)
    '''

    try:
        with pd.HDFStore(file, mode='r') as hdf5file:
            if "/" + key.strip("/") in hdf5file.keys():
                return True
            else:
                return False

    except Exception as exception:
        print(exception)

    return

def KeysOfStore(file):

    '''Get list of keys in HDF store file

    Parameters:
        file: str
            Path to HDF store file

    Returns
        list
           Keys
    
    Usage:
        KeysOfStore(file)
    '''

    try:
        with pd.HDFStore(file, mode='r') as hdf5file:
            return list(hdf5file.keys())

    except Exception as exception:
        print(exception)

    return

def downloadFile(url, saveDir, saveName = None):

    '''Download a file from internet

    Parameters:
        url: str
            URL to download from

        saveDir: str
            Path to save downloaded file to

        saveName: str, Default None
            New name for downloaded file

    Returns
        None
    
    Usage:
        downloadFile(url, saveDir)
    '''

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

    '''Read R data file

    Parameters:
        fullPath: str
            Path to file

        takeGeneSymbolOnly: boolean, Default True
            Wether to save gene symbol only

        saveToHDF: boolean, Default True
            Wether to save data in HDF format

        returnSizeOnly: boolean, Default False
            Get data size and return, without reading the data itself

    Returns
        None
    
    Usage:
        readRDataFile(path)
    '''

    if (not returnSizeOnly) and os.path.isfile(fullPath):
        print('File already exists:', fullPath)

        df = pd.read_hdf(fullPath, key='df')

        return df

    from rpy2.robjects import r as R

    R['load'](fullPath[:-len('.h5')])
    ls = np.array(R['ls']())

    rSparseMatrix = R[ls[0]]

    print('Matrix size:', end='\t', flush=True)
    size = R('dim')(rSparseMatrix)
    print(size)

    if returnSizeOnly:

        return size[0], size[1]

    columns = pd.Index(np.array(R['colnames'](rSparseMatrix))).astype(str)
    index = pd.Index(np.array(R['rownames'](rSparseMatrix))).astype(str)

    if takeGeneSymbolOnly:
        index = index.str.split('_ENS', expand=True).get_level_values(0)
    
    R['writeMM'](obj=rSparseMatrix, file='%s.mtx' % (fullPath))

    df = pd.DataFrame(index=index, columns=columns, data=scipy.io.mmread('%s.mtx' % (fullPath)).toarray().astype(int))

    os.remove('%s.mtx' % (fullPath))

    df = df.loc[~df.index.duplicated(keep='first')]
    df = df.T.loc[~df.T.index.duplicated(keep='first')].T

    if saveToHDF:
        df.to_hdf(fullPath, key='df', mode='a', complevel=4, complib='zlib')

        return df

    return df
      
def normSum1(data):

    w = np.nansum(data)
    if w != w or w == 0.:
        w = 1.

    return np.nan_to_num(data) / w

def silhouette(data, n_clusters, cluster_labels):

    u = np.unique(cluster_labels)
    n_clusters = len(u)

    df = pd.DataFrame(data=data, index=cluster_labels, columns=cluster_labels)
    print(df)

    print('DCS')
    dft = pd.DataFrame()
    for t in u:
        temp = df.loc[df.index==t, df.columns==t].values
        if temp.shape[0]>1:
            inside = temp[np.triu_indices(temp.shape[0], 1)].mean()
            outside = df.loc[df.index==t, df.columns!=t].values.mean()
            print(t, '\t', temp.shape[0], '\t', np.round(inside, 2), '\t', np.round(outside, 2))

    if type(cluster_labels[0]) in [str, np.str_]:
        dnames = {i:l for i,l in enumerate(u)}
        rdnames = {l:i for i,l in enumerate(u)}
        cluster_labels = np.array([rdnames[l] for l in cluster_labels])

    silhouette_avg = silhouette_score(data, cluster_labels)
    sample_silhouette_values = silhouette_samples(data, cluster_labels)

    fig, ax1 = plt.subplots(figsize=(7,7))

    y_lower = 1
    for i in range(n_clusters):
        ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]
        ith_cluster_silhouette_values.sort()
        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i
        color = cm.nipy_spectral(float(i) / n_clusters)
        ax1.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)
        try:
            label = str(dnames[i])
        except:
            label = str(i)

        ax1.text(0., y_lower + 0.5 * size_cluster_i, '' + str(dnames[i]) + ', avg=' + str(np.round(np.mean(ith_cluster_silhouette_values), 2)) + ', max=' + str(np.round(np.max(ith_cluster_silhouette_values), 2)), fontsize=10, color='k', va='center', ha='center').set_path_effects([path_effects.Stroke(linewidth=1., foreground='white'), path_effects.Normal()])
        y_lower = y_upper + 1

    ax1.set_xlabel("Silhouette coefficient values")
    ax1.set_xlim([min(sample_silhouette_values) - 0.1, max(sample_silhouette_values) + 0.1])
    ax1.set_ylim([0, len(data) + (n_clusters + 1) * 1])

    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")
    ax1.text(silhouette_avg + 0.025, ax1.get_ylim()[1], 'avg=' + str(np.round(silhouette_avg, 2)), color="red")

    ax1.spines['left'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.set_yticklabels([])
    ax1.set_yticks([])

    plt.show()

    return