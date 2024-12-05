import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import euclidean_distances

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype = {'subject_id': str, 'device_id': str})

distances = euclidean_distances(data[['x', 'y']], data[['x', 'y']])
distances[distances == 0] = np.nan

mean_per_row = np.nanmean(distances, axis=1)

eps = mean_per_row.std()
eps = pd.DataFrame({'eps': eps}, index=[0])

eps.to_csv(output_file, index=None, header=True)
