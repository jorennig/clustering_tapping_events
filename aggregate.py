import pandas as pd

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'subject_id': str})

cols = ['subject_id', 'test', 'feature_digital', 'two_week_period']
summary = data.groupby(cols)['value_digital'].mean().reset_index()

summary.to_csv(output_file, index=False)
