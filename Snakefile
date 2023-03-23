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
#DATASETS, SAMPLES, _ = glob_wildcards(os.path.join(OUTDIR, "{dataset,[^/]+}",
#        "{sample,[^/]+}", "{unused,[^/]+}.fasta"))
# run only successuful pharokka runs
DATASETS, SAMPLES = glob_wildcards(os.path.join(OUTDIR, "{dataset,[^/]+}",
        "{sample,[^/]+}", "pharokka", "phanotate.faa"))

if not os.path.isfile(os.path.join(OUTDIR, "extracted_ok")):
    ALL_OUTPUT = os.path.join(OUTDIR, "snakemake_genomes")
else:
    ALL_OUTPUT = expand(os.path.join(OUTDIR_SAMPLE, "defense_finder", "defense_finder_systems.tsv"), zip, dataset=DATASETS, sample=SAMPLES)

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
    threads:
        2
    resources:
        slurm_partition = "standard",
        runtime = 120,
        mem_mb = 8192,
        slurm_extra = "-J extract"
    params:
        filter="Provirus",
        length=1,
        completeness=100,
        outdir=OUTDIR
    shell:
        "scripts/extract_topology.py {params.filter} {params.length} {params.completeness} {input.tsv} {input.mfa} {params.outdir}"

rule snakemake_genomes:
    input:
        os.path.join(OUTDIR, "extracted_ok")
    output:
        os.path.join(OUTDIR, "snakemake_genomes")
    threads:
        2
    resources:
        slurm_partition = "long",
        runtime = 20160,
        mem_mb = 8192,
        slurm_extra = "-J genomes"
    shell:
        " ".join(sys.argv)

rule pharokka_annotation:
    input:
        os.path.join(OUTDIR_SAMPLE, "{sample}.fasta")
    output:
        os.path.join(OUTDIR_SAMPLE, "pharokka", "phanotate.faa")
    threads:
        32
    resources:
        slurm_partition = "standard",
        runtime = 7,
        mem_mb = 8192,
        slurm_extra = "-J pharokka  --exclusive=user"
    params:
        outdir=os.path.join(OUTDIR_SAMPLE, "pharokka")
    shell:
        "/bin/timeout -k 60s 30m pharokka.py -i {input} -o {params.outdir} -t {threads} -f"

rule defence_finder:
    input:
        os.path.join(OUTDIR_SAMPLE, "pharokka", "phanotate.faa")
    output:
        os.path.join(OUTDIR_SAMPLE, "defense_finder", "defense_finder_systems.tsv")
    threads:
        32
    resources:
        slurm_partition = "standard",
        runtime = 4,
        mem_mb = 8192,
        slurm_extra = "-J definder --exclusive=user"
    group:
        "pharokka"
    params:
        outdir=os.path.join(OUTDIR_SAMPLE, "defense_finder")
    log:
        os.path.join(OUTDIR_SAMPLE, "defense_finder_systems.log")
    shell:
        "/bin/timeout -k 60s 30m defense-finder run -o {params.outdir} {input} > {log}"
