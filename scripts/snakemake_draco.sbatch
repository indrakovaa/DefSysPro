#!/bin/bash
#SBATCH --job-name snake_draco
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8192
#SBATCH --output  log/snake_draco.%j.out               #%j means job ID
#SBATCH --error   log/snake_draci.%j.err
#SBATCH --time=7-00:00:00 #days
#SBATCH --exclusive=user                               # on node can be only user jobs

unset SLURM_MEM_PER_CPU SLURM_MEM_PER_GPU SLURM_MEM_PER_NODE

#modules
#module load tools/python/3.8

#Commands
cd ~/DefSystems
snakemake --unlock
snakemake extract_genomes -c all # runs Snakefile rule extract_genomes
#check if genome extraction was done, then decide 
if [ -f "output/extracted_ok" ]; then
        if [ "$1" == "local" ]; then
                snakemake -c all
        else
                # Quota = 10 nodes per user, we're using 3 jobs per node + one long
                # (-j) 27/3 = 9 (+ 1 long)
                snakemake --slurm -j 27 --max-status-checks-per-second 1 --keep-going \
                --group-components pharokka=100 definder=100 --resource cpus_per_task=32 mem_mb=8192
        fi
else
        echo "ERROR: No output/extracted_ok file exists = Missing correct genomes." >&2
fi
