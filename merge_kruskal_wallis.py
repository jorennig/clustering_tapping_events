import pandas as pd

input_file = snakemake.input.data[0]
input_kruskal_wallis = snakemake.input.kruskal_wallis[0]
cols = snakemake.params.group_cols
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'subject_id': str})
kruskal_wallis = pd.read_csv(input_kruskal_wallis, dtype={'subject_id': str})

cols = cols + ['n_clusters']
merged = data.merge(kruskal_wallis, on=cols, how='left')

merged.to_csv(output_file, index=False)
