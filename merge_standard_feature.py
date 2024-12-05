import pandas as pd

input_file = snakemake.input.data[0]
input_features = snakemake.input.features
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype = {'subject_id': str})
features = pd.read_csv(input_features, dtype = {'subject_id': str})

data = data[data['feature_name']=='intertap_interval_mean']

cols = ['', '']
data = data.merge(features, on=cols)

data.to_csv(output_file, index=None, header=True)
