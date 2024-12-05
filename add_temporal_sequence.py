import pandas as pd

input_file = snakemake.input.data[0]
cols = snakemake.params.group_cols
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype = {'subject_id': str, 'device_id': str})
data['sequence'] = ''

if not data.empty:
    data['sequence'] = data.groupby(cols).cumcount()

data.to_csv(output_file, index=None, header=True)
