import pandas as pd

traj_list = []
for file in snakemake.input:
    traj = pd.read_csv(file)
    traj_list.append(traj)

# combine into a single dataframe
agg_traj = pd.concat(traj_list, ignore_index=True)

# save output
agg_traj.to_csv(snakemake.output[0], index=False)