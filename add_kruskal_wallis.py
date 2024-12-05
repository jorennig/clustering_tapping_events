import pandas as pd

input_file = snakemake.input.data[0]
input_kruskal_wallis = snakemake.input.kruskal_wallis[0]
cols = snakemake.params.group_cols
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'subject_id': str})
kruskal_wallis = pd.read_csv(input_kruskal_wallis, dtype={'subject_id': str})

if (not data.empty) & (not kruskal_wallis.empty):
    cols = cols + ['n_clusters']
    data = data.merge(kruskal_wallis, on=cols, how='left') \
               .drop(columns=['stat', 'p'])

data.to_csv(output_file, index=False)
