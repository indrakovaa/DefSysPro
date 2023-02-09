#!/usr/bin/env python
import csv
import sys

from Bio import SeqIO

"""
  header fasta

  >IMGVR_UViG_2504643025_000001|2504643025|2504645545|381-5974

  name in the tsv file - informace v prvnim sloupci
  
  IMGVR_UViG_2504643025_000001

  python script to extract fasta entries matching criteria
"""

if (len(sys.argv) != 3):
  print("Usage: {} <tsv filename> <fasta filename>".format(sys.argv[0]))
  exit(1)

tsv = sys.argv[1]
fasta = sys.argv[2]

print("tsv: {}".format(tsv), file=sys.stderr)
print("fasta: {}".format(fasta), file=sys.stderr)

# Generate dictionary of seq_ids we're interested in
seq_ids = dict()
#with open(tsv, newline='') as csvfile:
with open(tsv) as csvfile:
    lines = csv.reader(csvfile, delimiter='\t')
    # Find UVIG and Topology index in tsv file header line
    header = next(lines)
    col_id = col_topology = -1
    for i, field in enumerate(header):
      if field == "UVIG":
        print("UVIG index: {}".format(i), file=sys.stderr)
        col_id = i
      elif field == "Topology":
        print("Topology index: {}".format(i), file=sys.stderr)
        col_topology = i
      else:
        continue
      if col_id >= 0 and col_topology >= 0:
        break
    # Process every tsv line
    for seq in lines:
      if seq[col_topology] == "GVMAG":
        seq_ids[seq[col_id]] = True

print("Start farsta lookup ...", file=sys.stderr)
# Go through fasta file and print seqs we're interested in
seq_all = SeqIO.parse(fasta, "fasta")
seq_gvmag = (seq for seq in seq_all if seq.id.split("|")[0] in seq_ids)
SeqIO.write(seq_gvmag, sys.stdout, "fasta")

print("Done.", file=sys.stderr)