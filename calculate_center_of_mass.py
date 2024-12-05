import pandas as pd

input_file = snakemake.input.data[0]
cols = snakemake.params.group_cols
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype = {'subject_id': str, 'device_id': str})

cols_c = ['subject_id', 'device_id', 'test', 'test_repeat', 'test_run_begin', 
          'test_run_type', 'Strokes Count', 'n_clusters', 'cluster_label', 'x', 'y']
center_of_mass = pd.DataFrame(columns=cols_c)

if not data.empty:
    data = data[~(data['cluster_label']=='noise')]
    
    data['n_clusters'] = data.groupby(cols)['cluster_label'].transform('nunique')
    data = data[data['n_clusters']>1]
    
    cols = cols + ['n_clusters', 'cluster_label']
    center_of_mass = data.groupby(cols)['x', 'y'].mean().reset_index()
    
center_of_mass.to_csv(output_file, index=None, header=True)
