import pandas as pd

input_file = snakemake.input.data[0]
device_id = str(snakemake.wildcards.device_id)
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype = {'subject_id': str, 'device_id': str})

data = data[data['device_id']==device_id]

data.to_csv(output_file, index=None, header=True)
