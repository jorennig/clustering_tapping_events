import pandas as pd

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

paths = pd.read_csv(input_file, dtype = {'subject_id': str, 'device_id': str, 'test_repeat': str})

cols = ['x', 'y', 'subject_id', 'device_id', 'test', 'test_repeat', 'test_run_begin', 
        'day_of_study', 'test_run_type', 
        'Strokes Count', 'intertap_distance_mean_pixel', 'ACTION_DOWN events']
data_tot = pd.DataFrame(columns=cols)

if not paths.empty:
    data_tot = pd.DataFrame()
    for index, row in paths.iterrows():
                
        data = pd.read_csv(row.folder_file, sep='\t', parse_dates=['timestamp', 'system_time'])
        
        data = data[data['type'].str.contains('DOWN')] \
                   .sort_values(by=['timestamp']) \
                   .filter(['x', 'y'])
        data['folder_file'] = row.folder_file
        data['ACTION_DOWN events'] = data.shape[0]
        
        data_tot = pd.concat([data_tot, data])

    data_tot = data_tot.merge(paths, on=['folder_file']) \
                       .filter(cols)

data_tot.to_csv(output_file, index=None, header=True)
