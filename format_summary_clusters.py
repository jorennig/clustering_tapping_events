import pandas as pd

input_file = snakemake.input.data[0]
features = snakemake.params.features
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype = {'subject_id': str, 'device_id': str})

id_vars = ['subject_id', 'device_id', 'test', 'test_repeat', 'test_run_begin', 
           'day_of_study', 'test_run_type']
value_vars = features
data = pd.melt(data, id_vars=id_vars, value_vars=value_vars,
               var_name='feature_digital', value_name='value_digital')

data['week'] = data['day_of_study'] // 7
data['two_week_period'] = data['day_of_study'] // 14

cols = ['subject_id', 'test', 'test_repeat', 'day_of_study']
data = data.sort_values(by=cols) \
           .drop_duplicates()

data.to_csv(output_file, index=None, header=True)
