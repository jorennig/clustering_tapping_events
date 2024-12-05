import pandas as pd

input_file = snakemake.input.data[0]
input_distance = snakemake.input.distance[0]
cols = snakemake.params.group_cols
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'subject_id': str})
distance = pd.read_csv(input_distance, dtype={'subject_id': str})

if (not data.empty) & (not distance.empty):
    cols_filter = cols + ['distance', 'n_clusters']
    distance_max = distance.groupby(cols).max().reset_index().filter(cols_filter)
    
    cols = cols + ['n_clusters']
    data = data.merge(distance_max, on=cols, how='left')

data.to_csv(output_file, index=False)
