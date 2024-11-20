from snakemake.utils import min_version
min_version("8.0.0")

configfile: "config.yml"

rule all:
    """
    Collect the main outputs of the workflow.
    """
    input:
        "results/tables/counts.tsv",
        "results/multiqc.html",
        "results/supplementary.html"

def get_sample_url(wildcards):
    return config["samples"][wildcards.sample_id]

rule get_SRA_by_accession:
    """
    Retrieve a single-read FASTQ file from a remote repository

    The fastq file is retrieved with wget and piped directly to the
    seqtk program which samples a number of reads defined by the
    max_reads parameter. The fastq output from seqtk is in turn piped
    to gzip and stored as a compressed *.fastq.gz file.

    The actual URL for the file is obtained from the config which
    requires that each sample_id is defined in the configfile as for
    example:

    sample_id: "https://url/to/file"
    """
    output:
        "data/{sample_id}.fastq.gz"
    log:
        "results/logs/get_SRA_by_accession/{sample_id}.log"
    params:
        max_reads = config["max_reads"],
        url = get_sample_url
    shell:
        """
        wget -q -o {log} -O - {params.url} | seqtk sample - {params.max_reads} | gzip -c > {output[0]}
        """

rule fastqc:
    """
    Run FastQC on a FASTQ file.
    """
    output:
        "results/{sample_id}_fastqc.html",
        "results/fastqc/{sample_id}_fastqc.zip"
    input:
        "data/{sample_id}.fastq.gz"
    shadow: "minimal"
    shell:
        """
        # Run fastQC and save the output to the current directory
        fastqc {input} -q -o .

        # Move the files which are used in the workflow
        mv {wildcards.sample_id}_fastqc.html {output[0]}
        mv {wildcards.sample_id}_fastqc.zip {output[1]}
        """

rule multiqc:
    """
    Aggregate all FastQC reports into a MultiQC report.
    """
    output:
        html="results/multiqc.html",
        stats="results/multiqc/multiqc_general_stats.txt"
    input:
        expand("results/fastqc/{sample_id}_fastqc.zip", sample_id = config["samples"].keys())
    log:
        "results/logs/multiqc/multiqc.log"
    shadow: "minimal"
    shell:
        """
        # Run multiQC and keep the html report
        multiqc -n multiqc.html {input} 2>{log}
        mv multiqc.html {output.html}
        mv multiqc_data/multiqc_general_stats.txt {output.stats}
        """

def get_fasta_path(wildcards):
    return config["genomes"][wildcards.genome_id]["fasta"]

rule get_genome_fasta:
    """
    Retrieve the sequence in fasta format for a genome.
    """
    output:
        "data/ref/{genome_id}.fa.gz"
    log:
        "results/logs/get_genome_fasta/{genome_id}.log"
    params:
        fasta_path = get_fasta_path
    shell:
        """
        wget {params.fasta_path} -O {output} -o {log}
        """

def get_gff_path(wildcards):
    return config["genomes"][wildcards.genome_id]["gff3"]

rule get_genome_gff3:
    """
    Retrieve annotation in gff3 format for a genome.
    """
    output:
        "data/ref/{genome_id}.gff3.gz"
    log:
        "results/logs/get_genome_gff3/{genome_id}.log"
    params:
        gff3_path = get_gff_path
    shell:
        """
        wget {params.gff3_path} -O {output} -o {log}
        """

rule index_genome:
    """
    Index a genome using Bowtie 2.
    """
    output:
        expand("results/bowtie2/{{genome_id}}.{n}.bt2",n=["1", "2", "3", "4","rev.1", "rev.2"])
    input:
        "data/ref/{genome_id}.fa.gz"
    log:
        "results/logs/index_genome/{genome_id}.log"
    shadow: "minimal"
    shell:
        """
        # Bowtie2 cannot use .gz, so unzip to a temporary file first
        gunzip -c {input} > tempfile
        bowtie2-build tempfile results/bowtie2/{wildcards.genome_id} >{log}
        """

rule align_to_genome:
    """
    Align a fastq file to a genome index using Bowtie 2.
    """
    output:
        # Here the sample_id wildcard is constrained with \w+ to match only
        # 'word characters', i.e. letters, numbers and underscore
        temp("results/bam/{sample_id,\\w+}.bam")
    input:
        "data/{sample_id}.fastq.gz",
        expand("results/bowtie2/{genome_id}.{substr}.bt2",
            genome_id=config["genome_id"],
            substr=["1", "2", "3", "4", "rev.1", "rev.2"])
    log:
        expand("results/logs/align_to_genome/{{sample_id}}_{genome_id}.log",
            genome_id = config["genome_id"])
    shell:
        """
        bowtie2 -x results/bowtie2/{config[genome_id]} -U {input[0]} > {output} 2>{log}
        """

rule sort_bam:
    """
    Sort a bam file.
    """
    output:
        "results/bam/{sample_id}.sorted.bam"
    input:
        "results/bam/{sample_id}.bam"
    shell:
        """
        samtools sort {input} > {output}
        """

rule generate_count_table:
    """
    Generate a count table using featureCounts.
    """
    output:
        "results/tables/counts.tsv",
        "results/tables/counts.tsv.summary"
    input:
        bams=expand("results/bam/{sample_id}.sorted.bam", sample_id = config["samples"].keys()),
        annotation=expand("data/ref/{genome_id}.gff3.gz", genome_id = config["genome_id"])
    log:
        "results/logs/generate_count_table.log"
    shell:
        """
        featureCounts -t gene -g gene_id -a {input.annotation} -o {output[0]} {input.bams} 2>{log}
        """

rule make_supplementary:
    """
    Generate Supplementary Material
    """
    output:
        "results/supplementary.html"
    input:
        counts = "results/tables/counts.tsv",
        summary_file = "results/tables/counts.tsv.summary",
        multiqc_file = "results/multiqc/multiqc_general_stats.txt",
    log:
        "results/logs/make_supplementary.log"
    params:
        SRR_IDs = config["sample_ids"],
        GSM_IDs = config["sample_ids_geo"]
    shell:
        """
        quarto render code/supplementary_material.qmd \
            -P counts_file:../{input.counts} \
            -P multiqc_file:../{input.multiqc_file} \
            -P summary_file:../{input.summary_file} \
            -P srr_ids:"{params.SRR_IDs}" \
            -P gsm_ids:"{params.GSM_IDs}"
        mv code/supplementary_material.html {output}
        """