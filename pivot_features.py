import pandas as pd

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype = {'subject_id': str, 'device_id': str})

cols = ['device_id', 'test_repeat', 'filename', 'extract_path', 'subject_id',
       'study', 'test', 'test_run_begin', 'test_run_type', 'day_of_study',
       'two_week_period', 'qc_pass', 'standard_feature']
data = pd.pivot_table(data, values='numeric_value', index=cols, columns='feature_name') \
         .reset_index()

data.to_csv(output_file, index=None, header=True)
