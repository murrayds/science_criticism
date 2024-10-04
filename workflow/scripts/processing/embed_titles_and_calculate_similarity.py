import json
import os

import pandas as pd
import numpy as np

from transformers import AutoTokenizer
from adapters import AutoAdapterModel
import torch

titles = pd.read_csv(snakemake.input[0])
titles = titles.dropna(subset=['citing_title'])

# Load the SPECTER2 model and tokenizer
tokenizer = AutoTokenizer.from_pretrained('allenai/specter2_base')
model = AutoAdapterModel.from_pretrained('allenai/specter2_base')

#load the adapter(s) as per the required task, provide an identifier for the adapter in load_as argument and activate it
model.load_adapter("allenai/specter2", source="hf", load_as="proximity", set_active=True)


from sklearn.metrics.pairwise import cosine_similarity


# Function to process and calculate cosine similarity for each group of CitedPaperId
def calculate_cosine_similarity_by_group(titles, batch_size=32):
    cosine_similarities = []
    
    for cited_paper_id, group in titles.groupby("CitedPaperId"):
        # Embed the cited_title only once
        cited_title = group["cited_title"].iloc[0] + tokenizer.sep_token
        cited_inputs = tokenizer(
            [cited_title],
            padding=True,
            truncation=True,
            return_tensors="pt",
            return_token_type_ids=False,
            max_length=512
        )
        cited_output = model(**cited_inputs)
        cited_embedding = cited_output.last_hidden_state[:, 0, :].detach().cpu().numpy()
        
        # Embed each citing_title
        citing_titles = group["citing_title"].tolist()
        citing_titles = [title + tokenizer.sep_token for title in citing_titles]
        
        for i in range(0, len(citing_titles), batch_size):
            batch = citing_titles[i:i + batch_size]
            inputs = tokenizer(
                batch,
                padding=True,
                truncation=True,
                return_tensors="pt",
                return_token_type_ids=False,
                max_length=512
            )
            output = model(**inputs)
            embeddings = output.last_hidden_state[:, 0, :].detach().cpu().numpy()
            
            # Calculate cosine similarity
            similarities = cosine_similarity(cited_embedding, embeddings)
            for j, citing_title in enumerate(batch):
                citing_paper_id = group["CitingPaperId"].iloc[i + j]
                cosine_similarities.append({
                    "cited_id": int(cited_paper_id),
                    "citing_id": int(citing_paper_id),
                    "cosine_similarity": float(similarities[0][j])
                })
    
    return cosine_similarities

# Calculate cosine similarities
cosine_similarities = calculate_cosine_similarity_by_group(titles)

# Convert list of dictionaries to DataFrame
cosine_similarity_df = pd.DataFrame(cosine_similarities)

# Save DataFrame to csv file
cosine_similarity_df.to_csv(snakemake.output[0], index=False)
