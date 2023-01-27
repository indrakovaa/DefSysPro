#!/bin/bash
#SBATCH --job-name snake_refseq
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --output  log/snake_refseq.out.%j               #%j means job ID
#SBATCH --error   log/snake_refseq.err.%j
#SBATCH --time=7-00:00:00 #days

unset SLURM_MEM_PER_CPU SLURM_MEM_PER_GPU SLURM_MEM_PER_NODE

#modules
#module load tools/python/3.8

#Commands
cd ~/DefSystems
snakemake --unlock
snakemake --slurm -j 24 --max-status-checks-per-second 1 --keep-going
