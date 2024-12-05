import pandas as pd

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype = {'subject_id': str, 'device_id': str}) \
         .drop(columns=['x', 'y']) \
         .drop_duplicates()

data.to_csv(output_file, index=None, header=True)