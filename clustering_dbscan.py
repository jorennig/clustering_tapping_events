import pandas as pd
from sklearn.cluster import DBSCAN

def dbscan(chunk, min_samples):    
    X = chunk[['x', 'y']]
    eps = chunk['intertap_distance_mean_pixel'].iloc[0]
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(X)
    labels = db.labels_
    chunk['cluster_label'] = labels
    chunk['cluster_label'] = (chunk['cluster_label'] + 1).astype(str)
    chunk['cluster_label'] = 'cluster_' + chunk['cluster_label']
    chunk['cluster_label'] = chunk['cluster_label'].replace('cluster_0', 'noise')
    return chunk

input_file = snakemake.input.data[0]
min_samples = snakemake.params.min_samples
cols = snakemake.params.group_cols
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype = {'subject_id': str, 'device_id': str})
data['cluster_label'] = ''

if not data.empty:
    data = data.groupby(cols).apply(dbscan, min_samples)
    
data.to_csv(output_file, index=None, header=True)
