import pickle
from collections import defaultdict
import re

# Define the dd() function for the defaultdict, which is needed to load the file...
def dd():
    return defaultdict(int)

dict = {}
for file in snakemake.input:
    # Load the first pickel file in the list
    with open(file, 'rb') as f:
        zscores = pickle.load(f)
    
    # Extract the year from the filename using regex
    year = int(re.search(r'journal_zscores_(\d{4})\.pickle', file).group(1))

    dict[year] = zscores

# Save the jcounts dictionary to a pickle file
with open(snakemake.output[0], 'wb') as f:
    pickle.dump(dict, f)