import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances

def distances_clusters(chunk):
    labels = list(chunk['cluster_label'].unique())
    distances = pd.DataFrame(euclidean_distances(chunk[['x', 'y']], chunk[['x', 'y']]), 
                             index=labels, columns=labels).unstack().reset_index(name='distance')
    distances['distance'] = round(distances['distance'], 2)
    distances = distances[~(distances['distance']==0)].drop_duplicates(subset='distance', keep='first')
    return distances

input_file = snakemake.input.data[0]
cols = snakemake.params.group_cols
ppi = snakemake.params.ppi
mm_per_inch = snakemake.params.mm_per_inch
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype = {'subject_id': str, 'device_id': str})
cols_c = ['subject_id', 'device_id', 'test', 'test_repeat', 'test_run_begin',
          'test_run_type', 'Strokes Count', 'level_0', 'level_1', 'distance']
distance_clusters = pd.DataFrame(columns=cols_c)

if not data.empty:
    data['x'] = data['x']*(mm_per_inch/ppi)
    data['y'] = data['y']*(mm_per_inch/ppi)
    
    cols = cols + ['n_clusters']
    cols_d = cols + ['distance']
    distance_clusters = data.groupby(cols).apply(distances_clusters).reset_index() \
                            .drop_duplicates(subset=cols_d) \
                            .drop(columns='level_8')

distance_clusters.to_csv(output_file, index=None, header=True)
