import pandas as pd

input_file = snakemake.input.data[0]
cols = snakemake.params.cols
cols_group = snakemake.params.cols_group
thresholds = snakemake.params.thresholds
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'subject_id': str})

# TODO: cut function for binning, cumsum for cumulative distribution
features = []
for t in thresholds:
    f = 'threshold=' + str(t) + '%'
    data[f] = data['percent_multistroke'] <= t
    features = features + [f]

data = pd.melt(data, id_vars=cols, value_vars=features,
                var_name='feature_name', value_name='numeric_value')

data['numeric_value'] = data['numeric_value'].replace({True: 'n_tests_above', 
                                                       False: 'n_tests_below'})

summary = data.groupby(cols_group + ['numeric_value'])['numeric_value'] \
              .count() \
              .reset_index(name='n')

summary = pd.pivot_table(summary, index=cols_group, values='n', 
                           columns='numeric_value') \
            .reset_index() \
            .fillna(0)

summary['percent_tests_with_multistroke_events'] = summary['n_tests_below'] / \
                     (summary['n_tests_below'] + summary['n_tests_above'])*100
summary = summary.rename(columns={'feature_name': 'threshold'})

summary.to_csv(output_file, index=None, header=True)
