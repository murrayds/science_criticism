import pandas as pd
import numpy as np
from itertools import combinations
from tqdm import tqdm
from collections import defaultdict
import pickle

df = pd.read_csv(snakemake.input[0])

# A simple function to mark the default value of the defaultdict
def dd():
    return defaultdict(int)


if snakemake.params.shuffle == True:
    grouped = df.groupby(['ReferenceYear'])
    
    # Function to shuffle ReferenceJournalId within each group
    def shuffle_group(group):
        group['ReferenceJournalId'] = np.random.permutation(group['ReferenceJournalId'].values)
        return group
    
    # Apply the shuffle function to each group
    shuffled_df = grouped[['CitingPaperId', 'ReferenceJournalId', "ReferenceYear"]].apply(shuffle_group)
    
    # Reset index to remove the grouping
    shuffled_df = shuffled_df.reset_index(drop=True)
    df = shuffled_df

jcounts = defaultdict(dd)

# Get a list containing institutions for each individual
journal_lists = df.groupby("CitingPaperId")["ReferenceJournalId"].apply(list)

# iterate through every "CitingPaperId"
# Group by CitingPaperId and iterate through each group
for sublist in tqdm(journal_lists, desc="Processing papers", unit="paper"):
    # Calculate combinations of JournalIds (allowing repeats)
    journal_pairs = list(combinations(sublist, 2))
    
    for pair in journal_pairs:
        jcounts[min(pair)][max(pair)] += 1


# Save the jcounts dictionary to a pickle file
with open(snakemake.output[0], 'wb') as f:
    pickle.dump(jcounts, f)