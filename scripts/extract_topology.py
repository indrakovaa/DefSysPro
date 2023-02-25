#!/usr/bin/env python
import csv
import sys

from Bio import SeqIO
from pathlib import Path

"""
  Python script to extract single fasta files based on information in
  extarnal tsv table. Fasta header matches UVIG.
  Topology of the viral sequence (Linear, Provirus, GVMAG...)
  Estimated completeness of the sequence (%)
  Length of the sequence (nt)
  test run ./extract_topology.py "Direct terminal repeat" 1 1 ../Phage_genomes/IMG_VR/IMGVR_all_Sequence_information-high_confidence.tsv ../Phage_genomes/IMG_VR/IMGVR_all_nucleotides-high_confidence.fna output
  """

if (len(sys.argv) != 7):
  print(f"Usage: {sys.argv[0]} <filtered Topology> <minimal Length> <minimal completeness [%] <tsv filename> <fasta filename> <output directory>")
  exit(1)

debug = False
#debug = True

filter_topology = sys.argv[1]
min_length = int(sys.argv[2])
min_completeness = float(sys.argv[3])
tsv = sys.argv[4]
fasta = sys.argv[5]
output_dir = sys.argv[6]

if debug:
  print(f"topology filter: {filter_topology}", file=sys.stderr)
  print(f"min length filter: {min_length}", file=sys.stderr)
  print(f"min completeness filter: {min_completeness}", file=sys.stderr)
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
    col_id = col_topology = col_completeness = col_length = -1
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
      elif field == "Length":
        if debug:
          print(f"Length index: {i}", file=sys.stderr)
        col_length = i
      else:
        continue
      if col_id >= 0 and col_topology >= 0 and \
        col_completeness >= 0 and col_length >= 0:
        
        break
    # Process every tsv line
    # 1. filter by topology
    # 2. filter by minimal length
    # 3. filter by completeness, drop NA values
    for seq in lines:
      if seq[col_topology] == filter_topology and \
        int(seq[col_length]) >= min_length and \
        not seq[col_completeness] == "NA" and \
        float(seq[col_completeness]) >= min_completeness:

        seq_ids[seq[col_id]] = True
        if debug:
          print(seq[col_id])

if debug:
  print(f"\nNr of found IDs: {len(seq_ids)}", file=sys.stderr)
  print("\nStart fasta lookup ...", file=sys.stderr)

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
