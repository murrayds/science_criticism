import pandas as pd

item_list = []
for file in snakemake.input:
    item = pd.read_csv(file)
    item_list.append(item)

# combine into a single dataframe
agg_items = pd.concat(item_list, ignore_index=True)

# save output
agg_items.to_csv(snakemake.output[0], index=False)