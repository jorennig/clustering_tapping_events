import pandas as pd

input_file = snakemake.input.data[0]
input_clinical = snakemake.input.clinical[0]
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'subject_id': str})
clinical = pd.read_csv(input_clinical, dtype={'subject_id': str})

data = data[data['two_week_period']==0]

data = data.merge(clinical, on=['subject_id']).drop_duplicates()

data.to_csv(output_file, index=None, header=True)
