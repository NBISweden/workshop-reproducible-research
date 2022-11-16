from snakemake.utils import min_version
min_version("5.3.0")

configfile: "config.yml"

rule all:
    """
    Collect the main outputs of the workflow.
    """
    input:
        "results/tables/counts.tsv",
        "results/multiqc.html",
        "results/rulegraph.png",
        "results/supplementary.html"

rule get_SRA_by_accession:
    """
    Retrieve a single-read FASTQ file from SRA (Sequence Read Archive) by run accession number.

    max_reads: Maximal number of reads to download for each sample.
    """
    output:
        "data/raw_internal/{sample_id}.fastq.gz"
    params:
        max_reads = config["max_reads"]
    version: "1.0"
    shell:
        """
        fastq-dump {wildcards.sample_id} -X {params.max_reads} --readids \
            --dumpbase --skip-technical --gzip -Z > {output}
        """

rule fastqc:
    """
    Run FastQC on a FASTQ file.
    """
    output:
        "results/{sample_id}_fastqc.html",
        "intermediate/{sample_id}_fastqc.zip"
    input:
        "data/raw_internal/{sample_id}.fastq.gz"
    version: "1.0"
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
        stats="intermediate/multiqc_general_stats.txt"
    input:
        expand("intermediate/{sample_id}_fastqc.zip", sample_id = config["sample_ids"])
    log:
        "results/logs/multiqc/multiqc.log"
    version: "1.0"
    shadow: "minimal"
    shell:
        """
        # Run multiQC and keep the html report
        multiqc -n multiqc.html {input} 2> {log}
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
        "data/raw_external/{genome_id}.fa.gz"
    log:
        "results/logs/get_genome_fasta/{genome_id}.log"
    params:
        fasta_path = get_fasta_path
    version: "1.0"
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
        "data/raw_external/{genome_id}.gff3.gz"
    log:
        "results/logs/get_genome_gff3/{genome_id}.log"
    params:
        gff3_path = get_gff_path
    version: "1.0"
    shell:
        """
        wget {params.gff3_path} -O {output} -o {log}
        """

rule index_genome:
    """
    Index a genome using Bowtie 2.
    """
    output:
        expand("intermediate/{{genome_id}}.{n}.bt2",n=["1", "2", "3", "4","rev.1", "rev.2"])
    input:
        "data/raw_external/{genome_id}.fa.gz"
    log:
        "results/logs/index_genome/{genome_id}.log"
    version: "1.0"
    shadow: "minimal"
    shell:
        """
        # Bowtie2 cannot use .gz, so unzip to a temporary file first
        gunzip -c {input} > tempfile
        bowtie2-build tempfile intermediate/{wildcards.genome_id} > {log}
        """

rule align_to_genome:
    """
    Align a fastq file to a genome index using Bowtie 2.
    """
    output:
        # Here the sample_id wildcard is constrained with \w+ to match only
        # 'word characters', i.e. letters, numbers and underscore
        temp("intermediate/{sample_id,\w+}.bam")
    input:
        fastq = "data/raw_internal/{sample_id}.fastq.gz",
        index = expand("intermediate/{genome_id}.{substr}.bt2",
            genome_id=config["genome_id"],
            substr=["1", "2", "3", "4", "rev.1", "rev.2"])
    log:
        expand("results/logs/align_to_genome/{{sample_id}}_{genome_id}.log",
            genome_id = config["genome_id"])
    version: "1.0"
    shell:
        """
        bowtie2 -x intermediate/{config[genome_id]} -U {input.fastq} > {output} 2>{log}
        """

rule sort_bam:
    """
    Sort a bam file.
    """
    output:
        "{prefix}.sorted.bam"
    input:
        "{prefix}.bam"
    version: "1.0"
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
        bams=expand("intermediate/{sample_id}.sorted.bam", sample_id = config["sample_ids"]),
        annotation=expand("data/raw_external/{genome_id}.gff3.gz", genome_id = config["genome_id"])
    log:
        "results/logs/generate_count_table.log"
    version: "1.0"
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
        multiqc_file = "intermediate/multiqc_general_stats.txt",
        rulegraph = "results/rulegraph.png"
    log:
        "results/logs/make_supplementary.log"
    params:
        SRR_IDs = config["sample_ids"],
        GSM_IDs = config["sample_ids_geo"],
        GSE_ID = config["series_id_geo"]
    shell:
        """
        echo 'rmarkdown::render("code/supplementary_material.Rmd", \
                  output_file="supplementary.html", \
                  output_dir="results/", \
                  params=list(counts_file="{input.counts}", \
                              multiqc_file="{input.multiqc_file}", \
                              summary_file="{input.summary_file}", \
                              rulegraph_file="{input.rulegraph}", \
                              SRR_IDs="{params.SRR_IDs}", \
                              GSM_IDs="{params.GSM_IDs}", \
                              GSE_ID="{params.GSE_ID}"))' \
        | R --vanilla > {log} 2>&1
        """

rule generate_rulegraph:
    """
    Generate a rulegraph for the workflow.
    """
    output:
        "results/rulegraph.png"
    shell:
        """
        snakemake --rulegraph --configfile config.yml | dot -Tpng > {output}
        """