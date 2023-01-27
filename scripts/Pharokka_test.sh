#!/bin/bash
#Pharokka commands to annotate phage genomes
PHIGENOME="/home/adela/DefSystems/Phage_genomes/\
5Jan2023_refseq_genomes.fa"
pharokka.py -i $PHIGENOME -o /home/adela/DefSystems/Pharokka_Refseq/ \
-d /home/adela/build/pharokka_databases/ -t 4 -f -m
# -i input genome fasta file
# -o output directory
# -f force, overwrites the output directory, use when outputdir exists
# -t threads
# -d  <path/to/database/>
# -m use for multiple contigs/metagonimcs data
