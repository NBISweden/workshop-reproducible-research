from snakemake.utils import min_version
min_version("8.0.0")

rule all:
    """
    Collect the main outputs of the workflow.
    """
    input:
        "results/tables/counts.tsv",
        "results/multiqc/multiqc.html"

def get_sample_url(wildcards):
    samples = {
        "SRR935090": "https://figshare.scilifelab.se/ndownloader/files/39539767",
        "SRR935091": "https://figshare.scilifelab.se/ndownloader/files/39539770",
        "SRR935092": "https://figshare.scilifelab.se/ndownloader/files/39539773"
    }
    return samples[wildcards.sample_id]

rule get_SRA_by_accession:
    """
    Retrieve a single-read FASTQ file
    """
    output:
        "data/{sample_id}.fastq.gz"
    params:
        url = get_sample_url
    shell:
        """
        curl -L -A "Mozilla/5.0" {params.url} | seqtk sample - 25000 | gzip -c > {output[0]}
        """

rule fastqc:
    """
    Run FastQC on a FASTQ file.
    """
    output:
        "results/fastqc/{sample_id}_fastqc.html",
        "results/fastqc/{sample_id}_fastqc.zip"
    input:
        "data/{sample_id}.fastq.gz"
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
        html = "results/multiqc/multiqc.html",
        stats = "results/multiqc/multiqc_general_stats.txt"
    input:
        "results/fastqc/SRR935090_fastqc.zip",
        "results/fastqc/SRR935091_fastqc.zip",
        "results/fastqc/SRR935092_fastqc.zip"
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
        "data/ref/NCTC8325.fa.gz"
    shell:
        """
        wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/dna//Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.toplevel.fa.gz -O {output}
        """

rule get_genome_gff3:
    """
    Retrieve annotation in gff3 format for a genome.
    """
    output:
        "data/ref/NCTC8325.gff3.gz"
    shell:
        """
        wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/gff3/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325//Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.37.gff3.gz -O {output}
        """

rule index_genome:
    """
    Index a genome using Bowtie 2.
    """
    output:
        "results/bowtie2/NCTC8325.1.bt2",
        "results/bowtie2/NCTC8325.2.bt2",
        "results/bowtie2/NCTC8325.3.bt2",
        "results/bowtie2/NCTC8325.4.bt2",
        "results/bowtie2/NCTC8325.rev.1.bt2",
        "results/bowtie2/NCTC8325.rev.2.bt2"
    input:
        "data/ref/NCTC8325.fa.gz"
    shell:
        """
        # Bowtie2 cannot use .gz, so unzip to a temporary file first
        gunzip -c {input} > tempfile
        bowtie2-build tempfile results/bowtie2/NCTC8325

        # Remove the temporary file
        rm tempfile
        """

rule align_to_genome:
    """
    Align a fastq file to a genome index using Bowtie 2.
    """
    output:
        "results/bam/{sample_id,\\w+}.bam"
    input:
        "data/{sample_id}.fastq.gz",
        "results/bowtie2/NCTC8325.1.bt2",
        "results/bowtie2/NCTC8325.2.bt2",
        "results/bowtie2/NCTC8325.3.bt2",
        "results/bowtie2/NCTC8325.4.bt2",
        "results/bowtie2/NCTC8325.rev.1.bt2",
        "results/bowtie2/NCTC8325.rev.2.bt2"
    shell:
        """
        bowtie2 -x results/bowtie2/NCTC8325 -U {input[0]} > {output}
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
        "results/tables/counts.tsv"
    input:
        bams = ["results/bam/SRR935090.sorted.bam",
                "results/bam/SRR935091.sorted.bam",
                "results/bam/SRR935092.sorted.bam"],
        annotation = "data/ref/NCTC8325.gff3.gz"
    shell:
        """
        featureCounts -t gene -g gene_id -a {input.annotation} -o {output} {input.bams}
        """