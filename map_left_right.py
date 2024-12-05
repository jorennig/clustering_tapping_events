import pandas as pd

input_file = snakemake.input.data[0]
demographics_file = snakemake.input.demographics
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype = {'subject_id': str})
demographics = pd.read_csv(demographics_file, dtype = {'subject_id': str}) \
                 .filter(['subject_id', 'handedness'])

data = data.merge(demographics, on=['subject_id'])

data['test_hand'] = 'RIGHT'
data.loc[(data['test']=='STT-NONDOMINANT') & (data['handedness']=='right'), 'test_hand'] = 'LEFT'
data.loc[(data['test']=='STT-DOMINANT') & (data['handedness']=='left'), 'test_hand'] = 'LEFT'

data.to_csv(output_file, index=None, header=True)
