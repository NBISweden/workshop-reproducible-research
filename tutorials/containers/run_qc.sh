# Make needed directories
mkdir -p data
mkdir -p results/fastqc

# Define URLs to fastq files for each sample
SRR935090="https://figshare.scilifelab.se/ndownloader/files/39539767"
SRR935091="https://figshare.scilifelab.se/ndownloader/files/39539770"
SRR935092="https://figshare.scilifelab.se/ndownloader/files/39539773"

# Download fastq files from remote repository and put in data directory
curl -L  "$SRR935090" -o data/SRR935090.fastq.gz
curl -L  "$SRR935091" -o data/SRR935091.fastq.gz
curl -L  "$SRR935092" -o data/SRR935092.fastq.gz

# Run fastqc and output to results directory
fastqc  data/SRR935090.fastq.gz --outdir=results/fastqc/
fastqc  data/SRR935091.fastq.gz --outdir=results/fastqc/
fastqc  data/SRR935092.fastq.gz --outdir=results/fastqc/
