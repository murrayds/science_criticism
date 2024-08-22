# Now, how do we get the Z-scores for particular papers? I think we probably need a separate file for our focus...???
import pickle 
from itertools import combinations

import pandas as pd
import numpy as np

from tqdm import tqdm

def dd():
    # Needed to loac the pickle file
    return defaultdict(int)

# Load the Z-scores
with open(snakemake.input.zscores[0], 'rb') as zf:
    z_scores = pickle.load(zf)

# Next, load the papers of interest
letters = pd.read_csv(snakemake.input.letters[0])[["id", "original_year"]]
letters = letters.rename(columns={"original_year": "year"})

articles = pd.read_csv(snakemake.input.nonletters[0])[["id", "year"]]

all_papers = pd.concat([letters, articles])


refs = pd.read_csv(snakemake.input.refs[0])
focus = refs.merge(
    all_papers,
    left_on='CitingPaperId',
    right_on='id',
    how='inner'
)[["CitingPaperId", "ReferenceJournalId", 'year']]


focus = focus[focus.year <= snakemake.params.year_max]
focus = focus[focus.year >= snakemake.params.year_min]

#
# Now, we will count all pairs of journals apearing in the reference list of
# these papers and identify their z-scores, with which we will calculate the 
# median and 10th percentile score. 
#
journal_lists = focus.groupby("CitingPaperId")

# Each paper is represented by a list of z-scores...
records = {}
# But we don't actually have to count them...we just need to extract the zscores
for _, grp in tqdm(journal_lists, desc="Processing papers", unit="paper", mininterval=5):
    # Calculate combinations of JournalIds (allowing repeats)
    sublist = list(grp.ReferenceJournalId)
    id = int(grp["CitingPaperId"].iloc[0])
    year = int(grp["year"].iloc[0])
    journal_pairs = list(combinations(sublist, 2))

    # get the z-score for each pair
    paper_zscore_list = []
    for pair in journal_pairs:
        paper_zscore_list.append(z_scores[year][min(pair)][max(pair)])

    paper_zscore_array = np.array(paper_zscore_list)
    if len(paper_zscore_list) > 0: # If a z-score was found...
        records[id] = {
            "Zscore_median": np.median(paper_zscore_array),
            "Zscore_10th": np.percentile(paper_zscore_array, 10),
            "pairs": len(journal_pairs)
        }

# Convert to a pandas dataframe
df_paper_zscores = pd.DataFrame.from_dict(records, orient="index").reset_index()
df_paper_zscores = df_paper_zscores[df_paper_zscores['Zscore_10th'] != np.inf]

# Rename the `index` column
df_paper_zscores = df_paper_zscores.rename(columns={"index": "id"})

# Save the output
df_paper_zscores.to_csv(snakemake.output[0], index=False)


