#!/usr/bin/env python
import os
from Bio import SeqIO
import numpy as np

# read list of directories with results
fasta_dirs =  []

# Open the file and read the content in a list ## funguje
with open('DF_phi_list.txt', 'r') as filehandle:
    directories = filehandle.readlines()
    for line in directories:
        # Remove linebreak which is the last character of the string
        fasta_dir = line[:-1]
        # Add item to the list
        fasta_dirs.append(fasta_dir)

print(fasta_dirs)

# print headers for each fasta_dir ## funguje
root_path = "/home/ca47yas/DefSystems/Phage_genomes/Refseq/GBid/"
headers = []
for item in fasta_dirs:
    path = os.path.join(root_path, item, item + ".fasta")
    with open(path, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            headers.append(record.description)

np.savetxt("Phages_with_DF.csv", 
           headers,
           delimiter =", ", 
           fmt ='% s')