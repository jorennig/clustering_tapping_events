import pandas as pd

input_file = snakemake.input.data
clinical_features = snakemake.params.clinical_features
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'subject_id': str})

data = data[(data['visit']=='baseline') & 
            (data['feature'].isin(clinical_features))]

data = data.rename(columns={'feature': 'feature_clinical', 
                            'value': 'value_clinical'}).drop(columns=['test'])

data.to_csv(output_file, index=None, header=True)
