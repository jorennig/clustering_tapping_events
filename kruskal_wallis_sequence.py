import pandas as pd
from scipy.stats import kruskal

def kruskal_wallis(chunk):
    labels = list(chunk['cluster_label'].unique())
    value_array = []
    for l in labels:
        cluster = chunk.loc[chunk['cluster_label']==l, 'sequence']
        value_array.append(cluster)
    s, p = kruskal(*value_array)
    return pd.Series([s, p], index=['stat', 'p'])

input_file = snakemake.input.data[0]
cols = snakemake.params.group_cols
p_threshold = snakemake.params.p_threshold
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'subject_id': str})
data = data[~(data['cluster_label']=='noise')]
data = data[data['n_clusters']>1]

kw_res = pd.DataFrame(columns=(cols + ['stat', 'p']))
if not data.empty:
    cols = cols + ['n_clusters']
    kw_res = data.groupby(cols).apply(kruskal_wallis).reset_index()
    kw_res['temporal_difference_clusters'] = (kw_res['p'] < p_threshold).astype(int)

kw_res.to_csv(output_file, index=None, header=True)
