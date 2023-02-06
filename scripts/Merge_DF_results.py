#!/usr/bin/env python
import pandas as pd
import os
import pathlib

# list of desired files
FILES = []
for root, dirs, files in os.walk("DefFinder_results/"):
    for file in files:
        if file.endswith('genes.tsv'):
            FILES.append(os.path.join(root, file))
            

# read all data frames
dfs = [pd.read_csv(file, sep='\t') for file in FILES]
# create one table
all_DF_res = pd.concat(dfs)


all_DF_res.to_csv('DF_result.tsv', sep='\t', header=True, index=False)
# Jeste bych potrebovala nazev slozky - NC_001416 pridat do prvniho sloupce