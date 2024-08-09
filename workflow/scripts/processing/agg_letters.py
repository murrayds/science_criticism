import pandas as pd

letters_list = []
for index in range(len(snakemake.input.letters)):
    # we want to merge these into a common dataframe for all journals.....
    letters = pd.read_csv(snakemake.input.letters[index])
    articles = pd.read_csv(snakemake.input.impacts[index])

    letters["RelativeTime"] = letters.letter_year - letters.original_year
    
    # filter to the desired timeframe
    letters = letters.query('original_year<2020')
    letters = letters.query('original_year>=2000')

    # filter letters pointing to retracted publications
    letters = letters[letters.retracted != True]

    # merge the letters dataframe with the paper_impact dataframe on the 'id' column
    letters_with_impact = pd.merge(
        letters, 
        articles.loc[:, articles.columns != 'venue'], # remove the redundant "venue" column
        left_on='original_id', 
        right_on="id"
    )
    
    # Drop unecessary columns
    letters_with_impact = letters_with_impact.drop(columns=["retracted", "year"])

    # For analyses centering on the original paper, we do not care whether they received 
    # multiple critical letters. Drop duplicates to simplify later work.
    letters_with_impact = letters_with_impact.drop_duplicates(subset='original_id')

    letters_list.append(letters_with_impact)

# combine them all into a single dataframe
agg_letters = pd.concat(letters_list, ignore_index=True)

# write the output
agg_letters.to_csv(snakemake.output[0], index=False)