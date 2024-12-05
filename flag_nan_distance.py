import pandas as pd

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'subject_id': str})

data = data.dropna(subset=['intertap_distance_mean_pixel'])

data.to_csv(output_file, index=False)
