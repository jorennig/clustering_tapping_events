import pandas as pd

# single script with unstack

input_file = snakemake.input.data[0]
threshold_strokes = int(snakemake.params.threshold_strokes)
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype = {'subject_id': str, 'device_id': str})

keep_values = data[(data['feature_name']=='Strokes Count') & 
                   (data['numeric_value']>threshold_strokes)] \
                   .drop(columns=['feature_name', 'numeric_value'])

cols = ['device_id', 'test_repeat', 'filename', 'extract_path', 'subject_id',
       'study', 'test', 'test_run_begin', 'day_of_study','two_week_period', 
       'qc_pass', 'standard_feature', 'handedness', 'test_hand']
data = data.merge(keep_values, on=cols)

data.to_csv(output_file, index=None, header=True)
