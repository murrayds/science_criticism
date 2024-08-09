import pandas as pd

fields_list = []
for file in snakemake.input:
    fields = pd.read_csv(file)
    fields_list.append(fields)

# combine into a single dataframe
agg_fields = pd.concat(fields_list, ignore_index=True)

# save output
agg_fields.to_csv(snakemake.output[0], index=False)