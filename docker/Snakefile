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
        "results/supplementary.pdf"

rule get_SRA_by_accession:
    """
    Retrieve a single-read FASTQ file from SRA (Sequence Read Archive) by run accession number.

    max_reads: Maximal number of reads to download for each sample.
    """
    output:
        "data/raw_internal/{sra_id}.fastq.gz"
    params:
        max_reads = config["max_reads"]
    version: "1.0"
    shell:
        """
        fastq-dump {wildcards.sra_id} -X {params.max_reads} --readids \
            --dumpbase --skip-technical --gzip -Z > {output}

        # This clears a cache where SRA Tools reserves a lot of space
        cache-mgr --clear >/dev/null 2>&1
        """

rule fastqc:
    """
    Run FastQC on a FASTQ file.
    """
    input:
        "data/raw_internal/{id}.fastq.gz"
    output:
        "results/{id}_fastqc.html",
        "intermediate/{id}_fastqc.zip"
    version: "1.0"
    shadow: "minimal"
    shell:
        """
        # Run fastQC and save the output to the current directory
        fastqc {input} -q -o .

        # Move the files which are used in the workflow
        mv {wildcards.id}_fastqc.html {output[0]}
        mv {wildcards.id}_fastqc.zip {output[1]}
        """

rule multiqc:
    """
    Aggregate all FastQC reports into a MultiQC report.
    """
    input:
        expand("intermediate/{id}_fastqc.zip", id = config["sample_ids"])
    output:
        html = "results/multiqc.html",
        stats = "intermediate/multiqc_general_stats.txt"
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

rule get_genome_fasta:
    """
    Retrieve the sequence in fasta format for a genome.
    """
    output:
        "data/raw_external/{genome_id}.fa.gz"
    params:
        fasta_path = lambda wildcards: config["genomes"][wildcards.genome_id]["fasta"]
    log:
        "results/logs/get_genome_fasta/{genome_id}.log"
    version: "1.0"
    shell:
        """
        wget {params.fasta_path} -O {output} -o {log}
        """

rule get_genome_gff3:
    """
    Retrieve annotation in gff3 format for a genome.
    """
    output:
        "data/raw_external/{genome_id}.gff3.gz"
    log:
        "results/logs/get_genome_gff3/{genome_id}.log"
    params:
        gff3_path = lambda wildcards: config["genomes"][wildcards.genome_id]["gff3"]
    version: "1.0"
    shell:
        """
        wget {params.gff3_path} -O {output} -o {log}
        """

rule index_genome:
    """
    Index a genome using Bowtie 2.
    """
    input:
        "data/raw_external/{genome_id}.fa.gz"
    output:
        expand("intermediate/{{genome_id}}.{n}.bt2", n = ["1","2","3","4"]),
        expand("intermediate/{{genome_id}}.rev.{n}.bt2", n = ["1","2"])
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
    input:
        fastq = "data/raw_internal/{sra_id}.fastq.gz",
        index = expand("intermediate/{genome_id}.{n}.bt2", genome_id = config["genome_id"], n = ["1","2","3","4"]),
        index_rev = expand("intermediate/{genome_id}.rev.{n}.bt2", genome_id = config["genome_id"], n = ["1","2"])
    output:
        temp("intermediate/{sra_id,\w+}.bam")
    log:
        expand("results/logs/align_to_genome/{{sra_id}}_{genome_id}.log", genome_id = config["genome_id"])
    version: "1.0"
    run:
        # This gives the base name for the genome index, i.e. "intermediate/some_id"
        # rather than "intermediate/some_id.*.bt2"
        indexBase = input.index[0].replace('.1.bt2','')
        shell("bowtie2 -x " + indexBase + " -U {input.fastq} > {output} 2> {log}")

rule sort_bam:
    """
    Sort a bam file.
    """
    input:
        "{prefix}.bam"
    output:
        "{prefix}.sorted.bam"
    version: "1.0"
    shell:
        """
        samtools sort {input} > {output}
        """

rule generate_count_table:
    """
    Generate a count table using htseq-count.
    """
    input:
        bams=expand("intermediate/{sra_id}.sorted.bam", sra_id = config["sample_ids"]),
        annotation=expand("data/raw_external/{genome_id}.gff3.gz", genome_id = config["genome_id"])
    output:
        "results/tables/counts.tsv"
    version: "1.0"
    shadow: "minimal"
    shell:
        """
        # htseq-count cannot use .gz, so unzip to a temporary file first
        gunzip -c {input.annotation} > tempfile

        # Save the count table as a temporary file and then prepend a header line
        # with the sample names
        htseq-count --format bam --type gene --additional-attr description --idattr gene_id {input.bams} tempfile > tempfile2
        echo '{input.bams}' | tr ' ' '\t' | cat - tempfile2 > {output}
        """

rule make_supplementary:
    """
    Generate Supplementary Material PDF
    """
    input:
        counts = "results/tables/counts.tsv",
        multiqc_file = "intermediate/multiqc_general_stats.txt",
        rulegraph = "results/rulegraph.png"
    output:
        "results/supplementary.pdf"
    params:
        SRR_IDs = config["sample_ids"],
        GSM_IDs = config["sample_ids_geo"]
    shell:
        """
        echo 'rmarkdown::render("code/supplementary_material.Rmd", \
                  output_file="supplementary.pdf", \
                  params=list(counts_file="{input.counts}", \
                              multiqc_file="{input.multiqc_file}", \
                              rulegraph_file="{input.rulegraph}", \
                              SRR_IDs="{params.SRR_IDs}", \
                              GSM_IDs="{params.GSM_IDs}"))' \
        | R --vanilla
        mv code/supplementary.pdf results/
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
