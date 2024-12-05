import pandas as pd

input_file = snakemake.input.data[0]
cols = snakemake.params.group_cols
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'subject_id': str})
add_cols = ['n_clusters', 'percent_noise']
data[add_cols] = ''

if not data.empty:
    
    data = data.drop(columns=add_cols)
    clusters = data[~(data['cluster_label']=='noise')]
    clusters = clusters.groupby(cols)['cluster_label'] \
                       .nunique() \
                       .reset_index(name='n_clusters')
    
    data = data.merge(clusters, on=cols, how='left')
    data['n_clusters'] = data['n_clusters'].fillna(0)
    
    data['noise'] = data['cluster_label']=='noise'
    data['percent_noise'] = data.groupby(cols)['noise'].transform('sum') / \
                            data.groupby(cols)['noise'].transform('count')
    data = data.drop(columns=['noise'])
    
data.to_csv(output_file, index=None, header=True)
