import pickle
from collections import defaultdict

# Define the dd() function for the defaultdict, which is needed to load the file...
def dd():
    return defaultdict(int)

# Load the first pickel file in the list
with open(snakemake.input[0], 'rb') as f:
    jcounts = pickle.load(f)

for next_file in snakemake.input[1:]:
    # Load the first pickel file in the list
    with open(next_file, 'rb') as f2:
        jcounts_next = pickle.load(f2)
    
    # Iterate through all combinations of keys in jcounts_next
    for key1 in jcounts_next:
        for key2 in jcounts_next[key1]:
            # Ensure key1 is always the smaller value
            
            # If the key doesn't exist in jcounts, add it
            if key1 not in jcounts:
                jcounts[key1] = defaultdict(int)
            
            # Add the value from jcounts_next to jcounts
            jcounts[key1][key2] += jcounts_next[key1][key2]

# Save the jcounts dictionary to a pickle file
with open(snakemake.output[0], 'wb') as f:
    pickle.dump(jcounts, f)