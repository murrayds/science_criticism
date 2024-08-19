import pandas as pd

fields_list = []
for file in snakemake.input:
    fields = pd.read_csv(file)
    fields_list.append(fields)

# combine into a single dataframe
agg_fields = pd.concat(fields_list, ignore_index=True)


field_hierarchy = pd.read_csv(snakemake.params[0])

# Merge field_hierarchy into agg_fields
agg_fields = pd.merge(agg_fields, field_hierarchy, left_on='field', right_on='ChildFieldOfStudyId', how='left')

# Rename FieldOfStudyId to field_level0
agg_fields = agg_fields.rename(columns={'FieldOfStudyId': 'field_level0'})

# Use coalesce to fill null values in field_level0 with the original field
agg_fields['field_level0'] = agg_fields['field_level0'].fillna(agg_fields['field'])

# Drop the unnecessary ChildFieldOfStudyId column if it exists
if 'ChildFieldOfStudyId' in agg_fields.columns:
    agg_fields = agg_fields.drop(columns=['ChildFieldOfStudyId'])

# save output
agg_fields.to_csv(snakemake.output[0], index=False)