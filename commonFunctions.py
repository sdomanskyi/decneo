from general import *
from genes import *

cleanListString = lambda c: str(list(c)).replace(' ', '').replace("'", '').replace(']', '').replace('[', '').replace(',', ', ')

def movingAverageCenteredLooped(a, halfWindowSize):

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
            a_avg = movingAverageCenteredLooped(a, n)
            print(np.sum(a_avg))
            plt.plot(np.arange(len(a)), a_avg)

        plt.show()
    '''

    n = halfWindowSize

    b = np.append(np.append(a[-n:], a), a[:n])
    s = np.cumsum(b) / (2.*n + 1.)
    s[2*n+1:] -= s[:-2*n-1]

    return s[2*n:]

def get_mean_std_cov_ofDataframe(df):

    se_mean = df.mean(axis=1)
    se_std = df.std(axis=1)

    df['1 mean'] = se_mean
    df['2 std'] = se_std
    df['3 cov'] = se_std / se_mean
    df = df[['1 mean', '2 std', '3 cov']]

    df.columns = ['mean', 'std', 'cov']

    df = df.sort_values('mean', axis=0, ascending=False)

    return df

def extractBootstrapVariability(variant, filePath = '', savePath = ''):

    df = pd.read_hdf(filePath, key='df').fillna(0).set_index('gene', append=True).droplevel('order')[variant].unstack('gene').fillna(0.).T

    writer = pd.ExcelWriter(savePath)
    for species in np.unique(df.columns.get_level_values('species')):
        get_mean_std_cov_ofDataframe(df.xs(species, level='species', axis=1)).to_excel(writer, species)

    writer.save()

    return

def getGenesOfPeak(se, heightCutoff = 0.5):

    se /= se.max()
    i0 = np.argmax(se.values)
    genes = [se.index[i0]]

    for i in range(i0 + 1, i0 + 2*10**3, 1):
        try:
            if se.iloc[i] >= heightCutoff:
                genes.append(se.index[i])
            else:
                break
        except:
            break

    for i in range(i0-1, 0, -1):
        try:
            if se.iloc[i] >= heightCutoff:
                genes.append(se.index[i])
            else:
                break
        except:
            break

    return genes

def get_df_distance(df, metric = 'correlation', genes = [], analyzeBy = 'batch', minSize = 10, groupBatches = True, pname = None, cutoff = 0.05):

    print('df_expr', '%1.1fMb'%(sys.getsizeof(df) / 1024**2), flush=True)

    genes = np.unique(genes)

    print('\tMetric:', metric, '\tAnalyzeBy:', analyzeBy, '\tminSize:', minSize)
    list_df_corr = []
    batches = []
    index = set()

    temp_batches = np.unique(df.columns.get_level_values(analyzeBy).values)[:]
    print('Number of batches:', len(temp_batches), flush=True)

    for i, batch in enumerate(temp_batches):
        df_sample = df.xs(key=batch, level=analyzeBy, axis=1, drop_level=False)
        df_sample = df_sample.loc[((df_sample > 0).sum(axis=1) / df_sample.shape[1]) >= cutoff]
        df_sample = df_sample.loc[df_sample.sum(axis=1) > 0.]

        if df_sample.shape[1] < minSize:

            continue

        print('\t', '%s-%s' % (i, batch), df_sample.shape, end=' ')

        temp_genes = df_sample.index.intersection(genes)
        print(len(temp_genes), end=' ', flush=True)

        measure = cdist(df_sample.loc[temp_genes].values, df_sample.values, metric=metric).T

        se_batch_corr = pd.DataFrame(measure, index=df_sample.index, 
                                     columns=temp_genes).stack().reorder_levels([1, 0]).sort_index()

        se_batch_corr = pd.concat([se_batch_corr], keys=[batch], axis=1, sort=False)
        se_batch_corr = pd.concat([se_batch_corr], keys=[metric], axis=1, sort=False)
        list_df_corr.append(se_batch_corr.copy())
        batches.append(batch)
        index = index.union(se_batch_corr.index.values)

    del df
    print()

    print('Making index', flush=True)
    index = pd.MultiIndex.from_tuples(list(index))
    print('Sorting index', flush=True)
    index = index.sort_values()
    print('\Index pairs:', len(index), '\tSelected genes:', len(np.unique(index.get_level_values(0))), '\tAll genes:', len(np.unique(index.get_level_values(1))))

    df_batches = pd.DataFrame(index=index, columns=batches, dtype=float)

    print('Re-indexing each sample', flush=True)
    for ibatch, batch in enumerate(batches):
        print('\t', ibatch, end=' ', flush=True)
        df_temp = list_df_corr[ibatch].replace(0., np.nan)
        df_batches.loc[df_temp.index, batch] = df_temp.values
        list_df_corr[ibatch] = None

    print()
    print(df_batches, flush=True)

    print('df_batches', '%1.1fMb'%(sys.getsizeof(df_batches) / 1024**2), flush=True)

    if groupBatches:
        df_batches = pd.Series(data=np.nanmedian(df_batches.values, axis=1), index=df_batches.index).unstack(0)

        for gene in df_batches.columns:
            try:
                df_batches.loc[gene, gene] = 0.
            except:
                pass

    return df_batches

def reduce(vIn, size = 100):

    v = vIn.copy()

    wh = np.where(~np.isnan(v))[0]

    if len(wh) > 0:
        bins =  np.linspace(np.min(v[wh]), np.max(v[wh]), num=size)

        v[wh] = bins[np.digitize(v[wh], bins) - 1]
    else:
        v[0] = 0.

    return v

def metric_euclidean_missing(u, v):

    wh = np.where(~np.isnan(u * v))[0]
                    
    return np.sqrt(((u[wh] - v[wh])**2).sum())

def nx_binom(nx_obj, enriched_genes, target_genes = False, background_genes = False):

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
            if not os.path.isfile('PCN.pklz'):
                nx_obj = nx.read_edgelist(nx_obj).to_undirected()
                write(nx_obj, 'PCN')
            else:
                nx_obj = read('PCN')
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
    
    data is np.array of shape (features, objects)

    method can be
    1 or 'Agglomerative'
    2 or 'Spectral'
    3 or 'KMeans'

    Methods 2 and 3 initial state is random, unless random_state is specified
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

    '''Test different number of clusters with KMeans, claclulate ARI for each'''

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
