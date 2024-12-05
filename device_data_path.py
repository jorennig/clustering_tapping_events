import pandas as pd
import os

input_file = snakemake.input.data[0]
device_data_path = snakemake.params.device_data_path
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype = {'subject_id': str, 'device_id': str, 
                                        'test_repeat': int})

data['folder_date'] = data['extract_path'].apply(lambda x: os.path.split(x)[1])

data['file_name'] = 'AT_' + data['device_id'] + '_' + data['folder_date'] + \
                    '_' + 'SPEED-TAPPING-' + data['test_repeat'].astype(str) + \
                    '_smartphone-touch.txt'

data['folder_file'] = data[['device_id', 'folder_date', 'file_name']] \
                          .apply(lambda row: os.path.join(device_data_path, *row), axis=1)

data = data.drop(columns=['extract_path', 'folder_date', 'file_name'])

data.to_csv(output_file, index=None, header=True)
