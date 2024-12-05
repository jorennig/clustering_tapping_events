import pandas as pd

input_file = snakemake.input.data[0]
min_proportion = snakemake.params.min_proportion
min_samples = snakemake.params.min_samples
cols = snakemake.params.group_cols
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'subject_id': str})

if not data.empty:
    
    cols = cols + ['cluster_label']
    data['n_datapoints'] = data.groupby(cols)['cluster_label'].transform('count')
    
    data.loc[data['n_datapoints'] < min_samples, 'cluster_label'] = 'noise'
    
    data['min_samples'] = round(data['ACTION_DOWN events']*min_proportion)
    data.loc[data['n_datapoints'] < data['min_samples'], 'cluster_label'] = 'noise'

data.to_csv(output_file, index=None, header=True)
