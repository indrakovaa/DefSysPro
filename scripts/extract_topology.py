#!/usr/bin/env python
import csv
import sys

from Bio import SeqIO
from pathlib import Path

"""
  header fasta

  >IMGVR_UViG_2504643025_000001|2504643025|2504645545|381-5974

  name in the tsv file - informace v prvnim sloupci
  
  IMGVR_UViG_2504643025_000001

  python script to extract fasta entries matching criteria
"""

if (len(sys.argv) != 5):
  print(f"Usage: {sys.argv[0]} <filtered Topology> <tsv filename> <fasta filename> <output directory>")
  exit(1)

#debug = False
debug = True

filter_topology = sys.argv[1]
tsv = sys.argv[2]
fasta = sys.argv[3]
output_dir = sys.argv[4]

if debug:
  print(f"topology filter: {filter_topology}", file=sys.stderr)
  print(f"tsv: {tsv}", file=sys.stderr)
  print(f"fasta: {fasta}", file=sys.stderr)
  print(f"output directory: {output_dir}\n", file=sys.stderr)

# Generate dictionary of seq_ids we're interested in
seq_ids = dict()
#with open(tsv, newline='') as csvfile:
with open(tsv) as csvfile:
    lines = csv.reader(csvfile, delimiter='\t')
    # Find UVIG and Topology index in tsv file header line
    header = next(lines)
    col_id = col_topology = col_completeness = -1
    for i, field in enumerate(header):
      if field == "UVIG":
        if debug:
          print(f"UVIG index: {i}", file=sys.stderr)
        col_id = i
      elif field == "Topology":
        if debug:
          print(f"Topology index: {i}", file=sys.stderr)
        col_topology = i
      elif field == "Estimated completeness":
        if debug:
          print(f"Completeness index: {i}", file=sys.stderr)
        col_completeness = i
      else:
        continue
      if col_id >= 0 and col_topology >= 0 and col_completeness >= 0:
        break
    # Process every tsv line
    for seq in lines:
      if seq[col_topology] == filter_topology and \
        (filter_topology in [ "GVMAG" ] or seq[col_completeness] == 100):
        
        seq_ids[seq[col_id]] = True

if debug:
  print("\nStart farsta lookup ...", file=sys.stderr)
# Go through fasta file and print seqs we're interested in
i = 0
for seq in SeqIO.parse(fasta, "fasta"):
  seq_id = seq.id.split("|")[0]
  if seq_id in seq_ids:
    j = i // 1000
    subdir = f"{j:06d}"
    i += 1
    filename = seq_id + ".fasta"
    if debug:
      print(f"j = {j:6d}, i = {i:6d}, filename: {filename}", file=sys.stderr)
    outdir_fasta = Path(output_dir, subdir, seq_id)
    outdir_fasta.mkdir(parents=True, exist_ok=True)
    SeqIO.write(seq, Path(outdir_fasta, filename), "fasta")

Path(output_dir, "extracted_ok").touch()
if debug:
  print(f"Done filtered {i} records.", file=sys.stderr)
