import pandas as pd
from glob import glob
import os

configfile: 'config.yaml'

dirs = sorted(glob(config['device_data_path'] + '/*'))
DEVICE_IDS = [os.path.basename(d) for d in dirs]
FEATURES = config['features']
FEATURES_CLUSTERING = config['features_clustering']
FEATURES_SUMMARY = config['features_summary']

rule select_features:
    input: script='scripts/select_features.py', 
           data=lambda w: config['feature_files'][w.test]
    output: 'data/{test}/select.csv'
    params: features=config['features']
    group: 'STT'
    script: 'scripts/select_features.py'
    
rule pivot_features:
    input: script='scripts/pivot_features.py',
           data=rules.select_features.output
    output: temp('data/{test}/pivot_features.csv')
    group: 'STT'
    script: 'scripts/pivot_features.py'

rule select_by_stroke_counts:
    input: script='scripts/select_by_threshold.py',
           data=rules.pivot_features.output
    output: 'data/{test}/select_by_stroke_counts.csv'
    params: feature='Strokes Count',
            threshold=config['threshold_strokes']
    group: 'STT'
    script: 'scripts/select_by_threshold.py'

rule flag_nan_distance:
    input: script='scripts/flag_nan_distance.py',
           data=rules.select_by_stroke_counts.output
    output: 'data/{test}/flag_nan_distance.csv'
    group: 'STT'
    script: 'scripts/flag_nan_distance.py'

rule device_data_path:
    input: script='scripts/device_data_path.py',
           data=rules.flag_nan_distance.output
    output: 'data/{test}/device_data_path.csv'
    params: device_data_path=config['device_data_path']
    group: 'STT'
    script: 'scripts/device_data_path.py'

rule concatenate:
    input: expand('data/{test}/device_data_path.csv', test=config['tests'])
    output: 'data/device_data_path.csv'
    run: pd.concat([pd.read_csv(f, dtype={'subject_id': str, 'device_id': str}) for f in input]).to_csv(output[0], index=False)

rule device_data_path_device_id:
    input: script='scripts/device_data_path_device_id.py',
           data=rules.concatenate.output
    output: 'data/device_id/{device_id}/device_data_path.csv'
    group: 'prepare_touch_data'
    script: 'scripts/device_data_path_device_id.py'

rule prepare_touch_data:
    input: script='scripts/prepare_touch_data.py',
           data=rules.device_data_path_device_id.output
    output: 'data/device_id/{device_id}/touch_data.csv'
    group: 'prepare_touch_data'
    script: 'scripts/prepare_touch_data.py'

rule clustering_dbscan:
    input: script='scripts/clustering_dbscan.py',
           data=rules.prepare_touch_data.output
    output: 'data/device_id/{device_id}/touch_data_clustered.csv'
    params: min_samples=config['min_samples'],
            group_cols=config['group_cols']
    group: 'clustering'
    script: 'scripts/clustering_dbscan.py'

rule flag_small_clusters:
    input: script='scripts/flag_small_clusters.py',
           data=rules.clustering_dbscan.output, 
    output: 'data/device_id/{device_id}/touch_data_n_samples.csv'
    params: min_proportion=config['min_proportion'],
            min_samples=config['min_samples'],
            group_cols=config['group_cols']
    group: 'clustering'
    script: 'scripts/flag_small_clusters.py'

rule n_clusters:
    input: script='scripts/n_clusters.py',
           data=rules.flag_small_clusters.output,
    output: 'data/device_id/{device_id}/touch_data_n_clusters.csv'
    params: group_cols=config['group_cols']
    group: 'clustering'
    script: 'scripts/n_clusters.py'

rule add_temporal_sequence:
    input: script='scripts/add_temporal_sequence.py',
           data=rules.n_clusters.output
    output: 'data/device_id/{device_id}/touch_data_sequence.csv'
    params: group_cols=config['group_cols']
    group: 'stats'
    script: 'scripts/add_temporal_sequence.py'

rule kruskal_wallis_sequence:
    input: script='scripts/kruskal_wallis_sequence.py',
           data=rules.add_temporal_sequence.output,
    output: 'data/device_id/{device_id}/kruskal_wallis_sequence.csv'
    params: group_cols=config['group_cols'],
            p_threshold=config['p_threshold']
    group: 'stats'
    script: 'scripts/kruskal_wallis_sequence.py'

rule concatenate_kruskal_wallis:
    input: expand('data/device_id/{device_id}/kruskal_wallis_sequence.csv', device_id=DEVICE_IDS)
    output: 'results/kruskal_wallis_clusters.csv'
    run: pd.concat([pd.read_csv(f, dtype={'subject_id': str, 'device_id': str}) for f in input]).sort_values(by=['subject_id', 'test', 'test_repeat', 'test_run_begin']).to_csv(output[0], index=False)

rule summary_clusters:
    input: script='scripts/summary_clusters.py',
           data=rules.n_clusters.output,
    output: 'data/device_id/{device_id}/summary_clusters.csv'
    group: 'summary'
    script: 'scripts/summary_clusters.py'

rule concatenate_summary_clusters:
    input: expand('data/device_id/{device_id}/summary_clusters.csv', device_id=DEVICE_IDS)
    output: ('data/summary_clusters_pre_format.csv')
    run: pd.concat([pd.read_csv(f, dtype={'subject_id': str, 'device_id': str}) for f in input]).to_csv(output[0], index=False)

rule merge_kruskal_wallis:
    input: script='scripts/merge_kruskal_wallis.py',
           data=rules.concatenate_summary_clusters.output,
           kruskal_wallis=rules.concatenate_kruskal_wallis.output,
    output: 'data/summary_kruskal_wallis.csv'
    params: group_cols=config['group_cols']
    group: 'summary'
    script: 'scripts/merge_kruskal_wallis.py'

rule format_summary_clusters:
    input: script='scripts/format_summary_clusters.py',
           data=rules.merge_kruskal_wallis.output
    output: 'data/summary_clusters.csv'
    params: features=config['features_clustering']
    group: 'summary'
    script: 'scripts/format_summary_clusters.py'

rule aggregate:
    input: script='scripts/aggregate.py',
           data=rules.format_summary_clusters.output
    output: 'data/summary_clusters_aggregated.csv'
    script: 'scripts/aggregate.py'

rule select_clinical:
    input: script='scripts/select_clinical.py',
           data=config['clinical']
    output: 'data/clinical_data.csv'
    params: clinical_features=config['clinical_features'],
    script: 'scripts/select_clinical.py'

rule merge_clinical:
    input: script='scripts/merge_clinical.py',
           data=rules.aggregate.output,
           clinical=rules.select_clinical.output,
    output: 'data/tapping_circle_clinical.csv'
    script: 'scripts/merge_clinical.py'

rule correlation_clinical:
    input: script='scripts/correlation.py',
           data=rules.merge_clinical.output
    output: 'results/correlation_clinical.csv'
    script: 'scripts/correlation.py'

rule scatterplot_clinical:
    input: script='scripts/scatterplot.py',
           data=rules.merge_clinical.output,
           results=rules.correlation_clinical.output
    output: 'results/scatterplot_clinical_baseline/{feature_clustering}.png'
    params: p_threshold=config['p_threshold_plot']
    script: 'scripts/scatterplot.py'

rule calculate_center_of_mass:
    input: script='scripts/calculate_center_of_mass.py',
           data=rules.flag_small_clusters.output
    output: 'data/device_id/{device_id}/center_of_mass_clusters.csv'
    params: group_cols=config['group_cols']
    script: 'scripts/calculate_center_of_mass.py'

rule distance_clusters:
    input: script='scripts/distance_clusters.py',
           data=rules.calculate_center_of_mass.output
    output: 'data/device_id/{device_id}/distance_clusters.csv'
    params: ppi=config['ppi'],
            mm_per_inch=config['mm_per_inch'],
            group_cols=config['group_cols']
    script: 'scripts/distance_clusters.py'

rule concatenate_distance_clusters:
    input: expand('data/device_id/{device_id}/distance_clusters.csv', device_id=DEVICE_IDS)
    output: 'data/distance_clusters_summary.csv'
    run: pd.concat([pd.read_csv(f, dtype={'subject_id': str, 'device_id': str}) for f in input]).to_csv(output[0], index=False)

rule merge_distance_clusters:
    input: script='scripts/merge_distance_clusters.py',
           data=rules.add_temporal_sequence.output,
           distance=rules.distance_clusters.output,
    output: 'data/device_id/{device_id}/touch_data_distance_clusters.csv'
    params: group_cols=config['group_cols']
    group: 'merge'
    script: 'scripts/merge_distance_clusters.py'

rule add_kruskal_wallis:
    input: script='scripts/add_kruskal_wallis.py',
           data=rules.merge_distance_clusters.output,
           kruskal_wallis=rules.kruskal_wallis_sequence.output,
    output: 'data/device_id/{device_id}/touch_data_kruskal_wallis.csv'
    params: group_cols=config['group_cols']
    group: 'merge'
    script: 'scripts/add_kruskal_wallis.py'

rule plot_clusters:
    input: script='scripts/plot_clusters.py',
           data=rules.add_kruskal_wallis.output,
           circle_coordinates=config['circle_coordinates'],
           icon_coordinates=config['icon_coordinates'],           
    output: touch('results/plots/junk/{device_id}/done')
    params: destination_folder_plots=config['destination_folder_plots'],
            group_cols=config['group_cols'],
            pixels_x=config['pixels_x'],
            pixels_y=config['pixels_y'],
            width=config['width'],
            height=config['height']
    group: 'plot_clusters'
    script: 'scripts/plot_clusters.py'

rule clusters_across_patient:
    input: script='scripts/summary_clustering.py',
           data=rules.merge_kruskal_wallis.output
    output: 'results/summary_clusters.csv',
    params: cols=config['group_cols'],
            cols_group=['test', 'feature_name']
    script: 'scripts/summary_clustering.py'

rule clusters_per_patient:
    input: script='scripts/summary_clustering.py',
           data=rules.merge_kruskal_wallis.output
    output: 'data/summary_clusters_patient.csv',
    params: cols=config['group_cols'],
            cols_group=['subject_id', 'test', 'feature_name']
    script: 'scripts/summary_clustering.py'

rule format_clusters_per_patient:
    input: script='scripts/format_clusters_per_patient.py',
           data=rules.clusters_per_patient.output
    output: 'data/summary_clusters_patient_format.csv',
    script: 'scripts/format_clusters_per_patient.py'

rule merge_summary_clinical:
    input: script='scripts/merge_clinical.py',
           data=rules.format_clusters_per_patient.output,
           clinical=rules.select_clinical.output
    output: 'data/summary_clusters_patient_clinical.csv'
    script: 'scripts/merge_clinical.py'
    
rule correlation_summary_clinical:
    input: script='scripts/correlation.py',
           data=rules.merge_summary_clinical.output
    output: 'results/correlation_summary_clinical.csv'
    script: 'scripts/correlation.py'

rule scatterplot_summary_clinical:
    input: script='scripts/scatterplot_summary.py',
           data=rules.merge_summary_clinical.output,
           results=rules.correlation_summary_clinical.output
    output: 'results/scatterplot_summary_clinical/{feature_summary}.png'
    params: p_threshold=config['p_threshold_plot']
    script: 'scripts/scatterplot_summary.py'


rule all:
    input: rules.concatenate.output,
           expand('data/device_id/{device_id}/device_data_path.csv', device_id=DEVICE_IDS),           
           rules.concatenate_kruskal_wallis.output,
           expand('data/device_id/{device_id}/summary_clusters.csv', device_id=DEVICE_IDS),
           rules.aggregate.output,
           expand('results/scatterplot_clinical_baseline/{feature_clustering}.png', feature_clustering=FEATURES_CLUSTERING),
           expand('data/device_id/{device_id}/distance_clusters.csv', device_id=DEVICE_IDS),
           rules.concatenate_distance_clusters.output,
           expand('data/device_id/{device_id}/touch_data_kruskal_wallis.csv', device_id=DEVICE_IDS),
           expand('results/plots/junk/{device_id}/done', device_id=DEVICE_IDS),
           rules.clusters_across_patient.output,
           expand('results/scatterplot_summary_clinical/{feature_summary}.png', feature_summary=FEATURES_SUMMARY),
           