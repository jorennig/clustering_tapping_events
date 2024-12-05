import pandas as pd

input_file = snakemake.input.data[0]
feature = snakemake.params.feature
threshold = int(snakemake.params.threshold)
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype = {'subject_id': str, 'device_id': str})

data = data[data[feature]>threshold]

data.to_csv(output_file, index=None, header=True)
