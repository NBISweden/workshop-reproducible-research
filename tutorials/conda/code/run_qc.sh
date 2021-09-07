# Make needed directories
mkdir -p data/raw_internal
mkdir -p intermediate/fastqc
mkdir -p results/fastqc

# Download fastq files using sra-tools and put in data/raw_internal/:
fastq-dump SRR935090 -X 12000 --gzip -Z > data/raw_internal/SRR935090.fastq.gz
fastq-dump SRR935091 -X 12000 --gzip -Z > data/raw_internal/SRR935091.fastq.gz
fastq-dump SRR935092 -X 12000 --gzip -Z > data/raw_internal/SRR935092.fastq.gz

# Run fastqc, put output zip in intermediate/fastq/ and html in results/fastqc:
fastqc  data/raw_internal/SRR935090.fastq.gz --outdir=intermediate/fastqc/
fastqc  data/raw_internal/SRR935091.fastq.gz --outdir=intermediate/fastqc/
fastqc  data/raw_internal/SRR935092.fastq.gz --outdir=intermediate/fastqc/
mv intermediate/fastqc/*html results/fastqc/
