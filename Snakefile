# Snakefile for the phage genome annotation and annotation 
# of defence systems using Defense Finder

# run Defense finder on gembase (more genomes - CDSs in one file)
defense-finder run -dbtype gembase esco_genomes.faa

rule all:
    input:
        expand("DefFinder_results/{genome}/defense_finder_systems.tsv",
        genome=GENOMES)

rule extract_genomes:
    input:
        tsv="Phage_genomes/IMG_VR/IMGVR_all_Sequence_information-high_confidence.tsv"
        mfa="Phage_genomes/IMG_VR/IMGVR_all_nucleotides-high_confidence.fna"
    output:
    shell:
        "module load tools/python/3.8; scripts/extract_entry_imgvr.py"
        " {input.tsv} {input.mfa}"

rule pharokka_annotation:
    input:
        "Phage_genomes/Refseq/GBid/{sample}/{sample}.fasta"   
    threads:
        32
    resources:
        slurm_partition = "short",
        runtime = 60,
        cpus_per_task = 32,
        tasks = 1,
        mem_mb = 8192,
        slurm_extra = "-J pharokka"
    params:
        outdir="Pharokka_results/{sample}/"
    output:
        "Pharokka_results/{sample}/phanotate.faa"
    shell:
        "module load tools/python/3.8; pharokka.py -i {input.fasta}"
        " -o {params.outdir} -t {threads} -f"

rule defence_finder:
    input:
        "Pharokka_results/{sample}/phanotate.faa"
    resources:
        slurm_partition = "short",
        runtime = 10,
        nodes = 8,
        cpus_per_task = 32,
        tasks = 1,
        mem_mb = 2048,
        slurm_extra = "-J definder"og:
    params:
        outdir="DefFinder_results/{sample}/"
    output:
        "DefFinder_results/{sample}/defense_finder_systems.tsv"
    log:
        log="DefFinder_results/{sample}/defense_finder_systems.log"
    shell:
        "module load tools/python/3.8; defense-finder run -o {params.outdir}"
        " {input} > {log}"
