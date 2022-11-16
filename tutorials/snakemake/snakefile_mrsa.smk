from snakemake.utils import min_version
min_version("5.3.0")

rule all:
    """
    Collect the main outputs of the workflow.
    """
    input:
        "results/tables/counts.tsv",
        "results/multiqc.html",
        "results/rulegraph.png"

rule get_SRA_by_accession:
    """
    Retrieve a single-read FASTQ file from SRA (Sequence Read Archive) by run accession number.
    """
    output:
        "data/raw_internal/{sample_id}.fastq.gz"
    shell:
        """
        fastq-dump {wildcards.sample_id} -X 25000 --readids \
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
        html = "results/multiqc.html",
        stats = "intermediate/multiqc_general_stats.txt"
    input:
        "intermediate/SRR935090_fastqc.zip",
        "intermediate/SRR935091_fastqc.zip",
        "intermediate/SRR935092_fastqc.zip"
    shell:
        """
        # Run multiQC and keep the html report
        multiqc -n multiqc.html {input}
        mv multiqc.html {output.html}
        mv multiqc_data/multiqc_general_stats.txt {output.stats}

        # Remove the other directory that multiQC creates
        rm -rf multiqc_data
        """

rule get_genome_fasta:
    """
    Retrieve the sequence in fasta format for a genome.
    """
    output:
        "data/raw_external/NCTC8325.fa.gz"
    shell:
        """
        wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/dna//Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.toplevel.fa.gz -O {output}
        """

rule get_genome_gff3:
    """
    Retrieve annotation in gff3 format for a genome.
    """
    output:
        "data/raw_external/NCTC8325.gff3.gz"
    shell:
        """
        wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/gff3/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325//Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.37.gff3.gz -O {output}
        """

rule index_genome:
    """
    Index a genome using Bowtie 2.
    """
    output:
        "intermediate/NCTC8325.1.bt2",
        "intermediate/NCTC8325.2.bt2",
        "intermediate/NCTC8325.3.bt2",
        "intermediate/NCTC8325.4.bt2",
        "intermediate/NCTC8325.rev.1.bt2",
        "intermediate/NCTC8325.rev.2.bt2"
    input:
        "data/raw_external/NCTC8325.fa.gz"
    shell:
        """
        # Bowtie2 cannot use .gz, so unzip to a temporary file first
        gunzip -c {input} > tempfile
        bowtie2-build tempfile intermediate/NCTC8325

        # Remove the temporary file
        rm tempfile
        """

rule align_to_genome:
    """
    Align a fastq file to a genome index using Bowtie 2.
    """
    output:
        "intermediate/{sample_id,\w+}.bam"
    input:
        "data/raw_internal/{sample_id}.fastq.gz",
        "intermediate/NCTC8325.1.bt2",
        "intermediate/NCTC8325.2.bt2",
        "intermediate/NCTC8325.3.bt2",
        "intermediate/NCTC8325.4.bt2",
        "intermediate/NCTC8325.rev.1.bt2",
        "intermediate/NCTC8325.rev.2.bt2"
    shell:
        """
        bowtie2 -x intermediate/NCTC8325 -U {input[0]} > {output}
        """

rule sort_bam:
    """
    Sort a bam file.
    """
    output:
        "intermediate/{sample_id}.sorted.bam"
    input:
        "intermediate/{sample_id}.bam"
    shell:
        """
        samtools sort {input} > {output}
        """

rule generate_count_table:
    """
    Generate a count table using featureCounts.
    """
    output:
        "results/tables/counts.tsv"
    input:
        bams = ["intermediate/SRR935090.sorted.bam", "intermediate/SRR935091.sorted.bam", "intermediate/SRR935092.sorted.bam"],
        annotation = "data/raw_external/NCTC8325.gff3.gz"
    shell:
        """
        featureCounts -t gene -g gene_id -a {input.annotation} -o {output} {input.bams}
        """

rule generate_rulegraph:
    """
    Generate a rulegraph for the workflow.
    """
    output:
        "results/rulegraph.png"
    shell:
        """
        snakemake --snakefile snakefile_mrsa.smk --config max_reads=0 --rulegraph | dot -Tpng > {output}
        """
