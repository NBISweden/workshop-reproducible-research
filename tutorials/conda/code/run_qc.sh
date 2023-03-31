# Make needed directories
mkdir -p data/raw_internal
mkdir -p intermediate/fastqc
mkdir -p results/fastqc

# Define URLs to fastq files for each sample
SRR935090="https://figshare.scilifelab.se/ndownloader/files/39539767"
SRR935091="https://figshare.scilifelab.se/ndownloader/files/39539770"
SRR935092="https://figshare.scilifelab.se/ndownloader/files/39539773"

# Download fastq files from remote repository and put in data/raw_internal/:
curl -L  "$SRR935090" -o data/raw_internal/SRR935090.fastq.gz
curl -L  "$SRR935091" -o data/raw_internal/SRR935091.fastq.gz
curl -L  "$SRR935092" -o data/raw_internal/SRR935092.fastq.gz

# Run fastqc, put output zip in intermediate/fastq/ and html in results/fastqc:
fastqc  data/raw_internal/SRR935090.fastq.gz --outdir=intermediate/fastqc/
fastqc  data/raw_internal/SRR935091.fastq.gz --outdir=intermediate/fastqc/
fastqc  data/raw_internal/SRR935092.fastq.gz --outdir=intermediate/fastqc/
mv intermediate/fastqc/*html results/fastqc/
