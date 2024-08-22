import pickle
import numpy as np
from collections import defaultdict

# Define the dd() function for the defaultdict, which is needed to load the file...
def dd():
    return defaultdict(int)

# Load the first pickel file in the list, this is the observed data
with open(snakemake.input[0], 'rb') as obs_f:
    jcounts_observed = pickle.load(obs_f)

# The remainder are the null networks...
jcounts_null = []
for f in snakemake.input[1:]:
    with open(f, 'rb') as null_f:
        jcounts_null.append(pickle.load(null_f))

# Now we calculate Z-scores...
z_scores = defaultdict(dd)
for key1 in jcounts_observed:
    for key2 in jcounts_observed[key1]:

        # get the list of null values
        null_values = []
        for i in range(len(snakemake.input) - 1): # for every iteration
            null_values.append(jcounts_null[i][key1][key2])

        with np.errstate(divide='ignore',invalid='ignore'):
            z_scores[key1][key2] = float(
                (jcounts_observed[key1][key2] - np.mean(null_values)) / np.std(null_values)
            )

# Save the jcounts dictionary to a pickle file
with open(snakemake.output[0], 'wb') as f:
    pickle.dump(z_scores, f)