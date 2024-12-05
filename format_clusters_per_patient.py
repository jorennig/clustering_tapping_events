import pandas as pd

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'subject_id': str})

data = data.rename(columns={'feature_name': 'feature_digital', 'percent': 'value_digital'})
data['feature_digital'] = 'percent_' + data['feature_digital']
data['two_week_period'] = 0

data.to_csv(output_file, index=None, header=True)
