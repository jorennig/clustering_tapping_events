import pandas as pd
from scipy.stats import spearmanr

def correlation_clinical(chunk):
    n = chunk.shape[0]
    r, p = spearmanr(chunk['value_digital'], 
                     chunk['value_clinical'], 
                     nan_policy='omit')
    return pd.Series([r, p, n], index=['rho', 'p', 'n'])

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype = {'subject_id': str})

corr = data.groupby(['test', 'feature_digital', 'feature_clinical']) \
           .apply(correlation_clinical) \
           .reset_index()

corr.to_csv(output_file, index=None, header=True)
