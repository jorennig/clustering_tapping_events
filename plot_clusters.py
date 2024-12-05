import pandas as pd
import plotnine as p9

input_file = snakemake.input.data[0]
input_circle_coordinates = snakemake.input.circle_coordinates
input_icon_coordinates = snakemake.input.icon_coordinates
cols = snakemake.params.group_cols
pixels_x = snakemake.params.pixels_x
pixels_y = snakemake.params.pixels_y
width = snakemake.params.width
height = snakemake.params.height
destination_folder_plots = snakemake.params.destination_folder_plots

data = pd.read_csv(input_file, dtype = {'subject_id': str, 'device_id': str, 
                                        'test_repeat': str, 'day_of_study': str})
data = data[data['n_clusters']>1]
circle_coordinates = pd.read_csv(input_circle_coordinates)
icon_coordinates = pd.read_csv(input_icon_coordinates)

if not data.empty:
    for chunk in data.groupby(cols):
        data_c = chunk[1]
                
        data_c['Strokes Count'] = data_c['Strokes Count'].astype(int).astype(str)
        
        title = 'subject_id: ' + data_c['subject_id'].iloc[0] + '\n' + \
                'distance_clusters: ' + str(round(data_c['distance'].iloc[0])) + ' mm' + '\n' + \
                'Temporal difference between clusters: ' + str(bool(data_c['temporal_difference_clusters'].iloc[0]))
                
        symbols = {'noise': '*', 'cluster_1': 'o', 'cluster_2': 's', 
                  'cluster_3': '^', 'cluster_4': 'v', 'cluster_5': 'D'}  
        
        p = p9.ggplot() \
           + p9.geom_point(circle_coordinates, p9.aes(x='x', y='y'), stroke=0, size=3) \
           + p9.geom_point(icon_coordinates, p9.aes(x='x_left', y='y'), stroke=0, size=2) \
           + p9.geom_point(icon_coordinates, p9.aes(x='x_right', y='y'), stroke=0, size=2) \
           + p9.geom_point(data_c, p9.aes(x='x', y='y', color='sequence', shape='cluster_label'), 
                           alpha=0.40, stroke=0, size=4) \
           + p9.xlim(0, pixels_x) \
           + p9.ggtitle(title) \
           + p9.scale_color_gradient(low='blue', high='red') \
           + p9.scale_shape_manual(values=symbols) \
           + p9.scale_y_reverse(limits=[pixels_y, 0])
        
        filename = '_'.join(list(data_c[cols].iloc[0]))
        
        output_file = destination_folder_plots + '/' + filename + '.png'
        p.save(output_file, dpi=300, width=width, height=height, verbose=False)
