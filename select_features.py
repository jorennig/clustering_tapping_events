import pandas as pd

input_file = snakemake.input.data
features = snakemake.params.features
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype = {'subject_id': str})

data = data[data['feature_name'].isin(features)]

data.to_csv(output_file, index=None, header=True)
