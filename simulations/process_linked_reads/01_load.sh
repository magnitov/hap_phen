# TELL-seq data for NA12878
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR10689418/SRR10689418
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR10689419/SRR10689419
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR10689420/SRR10689420
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR10689421/SRR10689421
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR10689422/SRR10689422

fastq-dump --split-files -F --gzip ./SRR10689418
fastq-dump --split-files -F --gzip ./SRR10689419
fastq-dump --split-files -F --gzip ./SRR10689420
fastq-dump --split-files -F --gzip ./SRR10689421
fastq-dump --split-files -F --gzip ./SRR10689422

zcat SRR10689418_1.fastq.gz SRR10689419_1.fastq.gz SRR10689420_1.fastq.gz SRR10689421_1.fastq.gz SRR10689422_1.fastq.gz > TELLseq_NA12878_R1.fastq
zcat SRR10689418_2.fastq.gz SRR10689419_2.fastq.gz SRR10689420_2.fastq.gz SRR10689421_2.fastq.gz SRR10689422_2.fastq.gz > TELLseq_NA12878_R2.fastq
zcat SRR10689418_3.fastq.gz SRR10689419_3.fastq.gz SRR10689420_3.fastq.gz SRR10689421_3.fastq.gz SRR10689422_3.fastq.gz > TELLseq_NA12878_I1.fastq
gzip TELLseq_NA12878_*.fastq

fastqc --threads 3 TELLseq_NA12878_*

# stLFR data for NA12878
wget -O stLFR_NA12878_R1.fastq.gz https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/stLFR/stLFR.split_read.1.fq.gz
wget -O stLFR_NA12878_R2.fastq.gz https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/stLFR/stLFR.split_read.2.fq.gz

fastqc --threads 2 stLFR_NA12878_*