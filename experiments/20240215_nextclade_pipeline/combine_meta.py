import pandas as pd
import os
os.chdir("/home/patrick/git/computationalPhylogenetics/experiments/20240215_nextclade_pipeline/data")

meta1 = pd.read_csv("deer_covid_meta.tsv", sep = "\t")
meta2 = pd.read_csv("omicron_meta.tsv", sep = "\t")

df = pd.concat([meta1, meta2], axis = 0 )#.to_csv("combined_meta.tsv", sep = "\t")
df = df.drop_duplicates()
df.to_csv("combined_meta.tsv", sep = "\t")

