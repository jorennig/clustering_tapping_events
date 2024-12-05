import pandas as pd

input_file = snakemake.input.data[0]
output_summary = snakemake.output.summary
output_summary_patient = snakemake.output.summary_patient

data = pd.read_csv(input_file, dtype = {'subject_id': str}) \
         .drop(columns=['cluster_label', 'n_datapoints']) \
         .drop_duplicates()
data['n_clusters>1'] = data['n_clusters'] > 1

cols = ['subject_id', 'device_id', 'test', 'test_repeat', 'test_run_begin', 
        'test_run_type', 'Strokes Count']
features = ['temporal_difference_clusters', 'n_clusters>1']
data = pd.melt(data, id_vars=cols, value_vars=features,
                var_name='feature_name', value_name='numeric_value') \
         .dropna()

data['numeric_value'] = data['numeric_value'].astype(bool).astype(str)

summary = data.groupby(['feature_name', 'test', 'numeric_value'])['numeric_value'] \
              .count() \
              .reset_index(name='n')

summary = pd.pivot_table(summary, index=['feature_name', 'test'], values='n', 
                           columns='numeric_value') \
            .reset_index()

summary['percent'] = summary['True']/(summary['True'] + summary['False'])*100



summary_patient = data.groupby(['subject_id', 'test', 'feature_name', 'numeric_value'])['subject_id'] \
                      .count() \
                      .reset_index(name='n')

summary_patient = pd.pivot_table(summary_patient, index=['subject_id', 'test', 'feature_name'], values='n', 
                                 columns='numeric_value') \
                    .reset_index() \
                    .fillna(0)

summary_patient['percent'] = summary_patient['True']/(summary_patient['True'] + \
                             summary_patient['False'])*100

summary_patient = summary_patient.filter(['subject_id', 'test', 'feature_name', 'percent']) \
                                 .rename(columns={'feature_name': 'feature_digital', 
                                                  'percent': 'value_digital'})
summary_patient['feature_digital'] = 'percent_' + summary_patient['feature_digital']

summary.to_csv(output_summary, index=None, header=True)
summary_patient.to_csv(output_summary_patient, index=None, header=True)
