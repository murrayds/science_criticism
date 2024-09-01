import pandas as pd

articles_list = []
for index in range(len(snakemake.input.letters)):
    # load the data
    letters = pd.read_csv(snakemake.input.letters[index])
    articles = pd.read_csv(snakemake.input.impacts[index])

    # Filter to those within the specified year
    articles = articles[articles.year < 2020]
    articles = articles[articles.year >= 2000]

    # Remove those articels with fewer than 5 references, they tend to be commentaries, 
    # news items, and other non research articles that should not be in the analysis
    articles = articles.query('ReferenceCount >= 5') # Filter Comment
    
    # Remove articles that have recieved a critical letter...
    articles = articles[~articles.id.isin(letters.original_id)]

    # Remove critical letters themselves
    articles = articles[~articles.id.isin(letters.letter_id)]

    articles_list.append(articles)
    
agg_articles = pd.concat(articles_list, ignore_index=True)

month = pd.read_csv(snakemake.input.month)
agg_articles = pd.merge(agg_articles, month, on="id", how = "left")


agg_articles.to_csv(snakemake.output[0], index=False)