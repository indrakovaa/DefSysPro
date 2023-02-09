#!/usr/bin/env python
from Bio import SeqIO
from pathlib import Path

multifasta = "/home/ca47yas/DefSystems/Phage_genomes/Refseq/5Jan2023_refseq_genomes.fa"
outdir = "/home/ca47yas//DefSystems/Phage_genomes/Refseq/single_refseq_genomes/GBid/"

#Path(outdir).mkdir(parents=True, exist_ok=True)
for seq_record in SeqIO.parse(multifasta, "fasta"):
     filename = seq_record.id+".fasta"
     outdir_fasta = outdir+seq_record.id+"/"
     Path(outdir_fasta).mkdir(parents=True, exist_ok=True)
     #print(outdir_fasta+filename)
     SeqIO.write(seq_record, outdir_fasta+filename, "fasta")

## For large multifasta might be necessary to use indexing
 #record_dict = SeqIO.index("Fasta/f002", "fasta")
 #len(record_dict) # length of the dictionary
 #print(len(record_dict["gi|1348917|gb|G26685|G26685"])) #length of the record
#record_dict.close()

