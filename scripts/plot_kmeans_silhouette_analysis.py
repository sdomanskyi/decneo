from scRegulation.commonFunctions import *

if __name__=='__main__':
    
    from sklearn.datasets import make_blobs
    data = make_blobs(n_samples=500, n_features=2, centers=4, cluster_std=1, center_box=(-10.0, 10.0), shuffle=True, random_state=1)[0]
    n_clusters = 4
    cluster_labels = KMeans(n_clusters=n_clusters, random_state=10).fit_predict(data)

    silhouette(data, cluster_labels, n_clusters)