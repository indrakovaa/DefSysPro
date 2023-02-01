#!/usr/bin/env python
import pandas as pd
import os

all_DF_res = pd.DataFrame()

# open desired file, append to new DF
for root, dirs, files in os.walk("DefFinder_results/"):
    for file in files:
        if file.endswith('genes.tsv'):
            data = pd.read_csv(file, sep='\t')
            all_DF_res = all_DF_res.append(data)

all_DF_res.to_csv('DF_result.tsv', sep='\t', index=False)