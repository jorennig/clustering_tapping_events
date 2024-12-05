import pandas as pd
from plotnine import ggplot, aes, geom_point, facet_wrap, ylab

input_file = snakemake.input.data[0]
input_results = snakemake.input.results[0]
feature = snakemake.wildcards.feature_clustering
p_threshold = snakemake.params.p_threshold
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'subject_id': str})
data = data[data['feature_digital']==feature]

results = pd.read_csv(input_results)

results['rho_str'] = 'rho=' + results['rho'].round(2).astype(str)
results['p_str'] = 'p=' + results['p'].round(4).astype(str)
results.loc[results['p']<p_threshold, 'p_str'] = 'p<' + str(p_threshold)
results['header'] = results['test'] + '\n ' + results['feature_clinical'] \
                    + '\n ' + results['rho_str'] + ', ' + results['p_str']

data = pd.merge(data, results, on=['test', 'feature_clinical', 'feature_digital'])

width = 15
height = 5
plot = ggplot(data, aes(x='value_clinical', y='value_digital')) \
        + geom_point(alpha=0.30, stroke=0, size=3) \
        + facet_wrap('header', nrow=2, scales='free_x') \
        + ylab(feature)

plot.save(output_file, dpi=300, width=width, height=height, verbose=False)
