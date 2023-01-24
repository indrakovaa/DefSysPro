# Snakefile for the phage genome annotation and annotation 
# of defence systems using Defense Finder and Padloc

# manually download genomes and create genome list
with open("../Phage_genomes/"
    "Refseq/single_refseq_genomes/40_GBid/10genome_names.txt") as f:
 GENOMES = f.read().splitlines()

# input configfile with all phage genomes
configfile: "config10.yaml"
# define function to define inputs in later stage of execution of Snakefile
def get_phage_genome_fasta(wildcards): # takes wildcards object
    return config["genomes"][wildcards.sample] #access it via attributes

rule all:
 input:
  expand("../DefFinder_results/{genome}/defense_finder_systems.tsv",
   genome=GENOMES)

rule pharokka_annotation:
    input:
        dbdir="../build/pharokka_databases/",
        fasta=get_phage_genome_fasta 
    threads:
        4
    params:
        outdir="../Pharokka_results/{sample}/"
    log:
        "../Pharokka_results/{sample}/{sample}.phanotate.log"
    output:
        "../Pharokka_results/{sample}/phanotate.faa"
    shell:
        "pharokka.py -i {input.fasta} -o {params.outdir}"
        " -d {input.dbdir} -t {threads} -f > {log}"

rule defence_finder:
    input:
        "../Pharokka_results/{sample}/phanotate.faa"
    log:
        log="../DefFinder_results/{sample}/defense_finder_systems.log"
    params:
        outdir="../DefFinder_results/{sample}/"
    output:
        "../DefFinder_results/{sample}/defense_finder_systems.tsv"
    shell:
        "defense-finder run -o {params.outdir} {input} > {log}"
