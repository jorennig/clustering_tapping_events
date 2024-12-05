import pandas as pd

input_file = snakemake.input.data[0]
cols = snakemake.params.cols
cols_group = snakemake.params.cols_group
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype = {'subject_id': str}) \
         .drop(columns=['cluster_label', 'n_datapoints']) \
         .drop_duplicates()
data['n_clusters>1'] = data['n_clusters'] > 1

features = ['temporal_difference_clusters', 'n_clusters>1']
data = pd.melt(data, id_vars=cols, value_vars=features,
                var_name='feature_name', value_name='numeric_value') \
         .dropna()

data['numeric_value'] = data['numeric_value'].astype(bool).astype(str)

summary = data.groupby(cols_group + ['numeric_value'])['numeric_value'] \
              .count() \
              .reset_index(name='n')

summary = pd.pivot_table(summary, index=cols_group, values='n', 
                           columns='numeric_value') \
            .reset_index() \
            .fillna(0)

summary['percent'] = summary['True']/(summary['True'] + summary['False'])*100

summary.to_csv(output_file, index=None, header=True)
