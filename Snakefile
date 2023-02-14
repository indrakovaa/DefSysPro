# Snakefile for the phage genome annotation and annotation 
# of defence systems using Defense Finder
import os.path
import sys

# Define output directories
OUTDIR = "output"
OUTDIR_SAMPLE = os.path.join(OUTDIR, "{dataset}/{sample}")

# Automatically detect what should all do
# either generate list of genomes to process
# or process the genomes
DATASETS, SAMPLES, _ = glob_wildcards(OUTDIR_SAMPLE + "/{unused}.fasta")
if not os.path.isfile(os.path.join(OUTDIR, "extracted_ok")):
    ALL_OUTPUT = os.path.join(OUTDIR, "snakemake_genomes")
else:
    ALL_OUTPUT = expand(os.path.join(OUTDIR_SAMPLE, "defense_finder_systems.tsv"), zip, dataset=DATASETS, sample=SAMPLES)

# Prepend loading correct python in every slurm host
if ("--slurm" in sys.argv):
    shell.prefix("module load tools/python/3.8; ")
else:
    shell.prefix("echo 'Add on slurm: module load tools/python/3.8;'; ")

rule all:
    input:
        ALL_OUTPUT

rule extract_genomes:
    input:
        tsv="Phage_genomes/IMG_VR/IMGVR_all_Sequence_information-high_confidence.tsv",
        mfa="Phage_genomes/IMG_VR/IMGVR_all_nucleotides-high_confidence.fna"
    output:
        os.path.join(OUTDIR, "extracted_ok")
    params:
        filter="GVMAG",
        outdir=OUTDIR
    shell:
        "scripts/extract_topology.py {params.filter} {input.tsv} {input.mfa} {params.outdir}"

rule snakemake_genomes:
    input:
        os.path.join(OUTDIR, "extracted_ok")
    output:
        os.path.join(OUTDIR, "snakemake_genomes")
    shell:
        " ".join(sys.argv)

rule pharokka_annotation:
    input:
        os.path.join(OUTDIR_SAMPLE, "{sample}.fasta")
    output:
        os.path.join(OUTDIR_SAMPLE, "phanotate.faa")
    threads:
        32
    resources:
        slurm_partition = "short",
        runtime = 60,
        mem_mb = 8192,
        slurm_extra = "-J pharokka"
    params:
        outdir=OUTDIR_SAMPLE
    shell:
        "pharokka.py -i {input} -o {params.outdir} -t {threads} -f"

rule defence_finder:
    input:
        os.path.join(OUTDIR_SAMPLE, "phanotate.faa")
    output:
        os.path.join(OUTDIR_SAMPLE, "defense_finder_systems.tsv")
    threads:
        32
    resources:
        slurm_partition = "short",
        runtime = 10,
        mem_mb = 2048,
        slurm_extra = "-J definder"
    params:
        outdir=OUTDIR_SAMPLE
    log:
        os.path.join(OUTDIR_SAMPLE, "defense_finder_systems.log")
    shell:
        "defense-finder run -o {params.outdir} {input} > {log}"
