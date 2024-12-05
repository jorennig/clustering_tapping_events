import pandas as pd
import numpy as np

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

input_file = r'\\exports.hps.kau.science.roche.com\pred\rpmda\BN40423\device_data\live\extracts\test-runs\353698103030148/2020-03-09_01-11-14/AT_353698103030148_2020-03-09_01-11-14_SPEED-TAPPING-1_smartphone-touch.txt'

input_file = r'\\exports.hps.kau.science.roche.com\pred\rpmda\BN40423\device_data\live\extracts\test-runs\353698103030148/2020-03-09_01-11-14/AT_353698103030148_2020-03-09_01-11-14_SPEED-TAPPING-1_smartphone-touch.txt'
output_file = r'\\exports.hps.kau.science.roche.com\pred\rpmda\BN40423-v2\jira\RPMDA-9719-visualize-raw-tapping-events\data\touch_data.csv'

data = pd.read_csv(input_file, sep='\t', parse_dates=['timestamp', 'system_time'])

n_events = len(data[data['type']=='ACTION_DOWN'])
seq_events = list(range(1, n_events+1))

data.loc[data['type']=='ACTION_DOWN', 'seq'] = seq_events
data['seq'] = data['seq'].ffill()

data['even_odd'] = data['seq'].mod(2).eq(0)

data['seconds'] = (data['system_time'] - data['system_time']).dt.total_seconds()

data.to_csv(output_file, index=None, header=True)
