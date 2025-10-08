# NA12878 LCLs from Genome in a Bottle:
wget -O NA12878_phased_possorted_bam.bam ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/10Xgenomics_ChromiumGenome_LongRanger2.1_09302016/NA12878_hg19/NA12878_hg19_phased_possorted_bam.bam
/DATA/users/m.magnitov/software/bamtofastq_linux --nthreads 16 NA12878_phased_possorted_bam.bam /DATA/users/m.magnitov/hap_phen/10X/fastq/NA12878/

# Yoruba Trio LCLs from 10X (NA19238, NA19239, NA19240):
wget -O NA19238_phased_possorted_bam.bam https://sra-pub-src-1.s3.amazonaws.com/ERR1868756/NA19238_phased_possorted_bam.bam.1
wget -O NA19239_phased_possorted_bam.bam https://sra-pub-src-1.s3.amazonaws.com/ERR1868757/NA19239_phased_possorted_bam.bam.1
wget -O NA19240_phased_possorted_bam.bam https://sra-pub-src-1.s3.amazonaws.com/ERR1868758/NA19240_phased_possorted_bam.bam.1
/DATA/users/m.magnitov/software/bamtofastq_linux --nthreads 16 NA19238_phased_possorted_bam.bam /DATA/users/m.magnitov/hap_phen/10X/fastq/NA19238/
/DATA/users/m.magnitov/software/bamtofastq_linux --nthreads 16 NA19239_phased_possorted_bam.bam /DATA/users/m.magnitov/hap_phen/10X/fastq/NA19239/
/DATA/users/m.magnitov/software/bamtofastq_linux --nthreads 16 NA19240_phased_possorted_bam.bam /DATA/users/m.magnitov/hap_phen/10X/fastq/NA19240/

# Han Chinese Trio LCLs from 10X (HG00512, HG00513, HG00514):
wget -O HG00512_phased_possorted_bam.bam https://sra-pub-src-1.s3.amazonaws.com/ERR1868750/HG00512_phased_possorted_bam.bam.1
wget -O HG00513_phased_possorted_bam.bam https://sra-pub-src-1.s3.amazonaws.com/ERR1868751/HG00513_phased_possorted_bam.bam.1
wget -O HG00514_phased_possorted_bam.bam https://sra-pub-src-1.s3.amazonaws.com/ERR1868752/HG00514_phased_possorted_bam.bam.1
/DATA/users/m.magnitov/software/bamtofastq_linux --nthreads 16 HG00512_phased_possorted_bam.bam /DATA/users/m.magnitov/hap_phen/10X/fastq/HG00512
/DATA/users/m.magnitov/software/bamtofastq_linux --nthreads 16 HG00513_phased_possorted_bam.bam /DATA/users/m.magnitov/hap_phen/10X/fastq/HG00513
/DATA/users/m.magnitov/software/bamtofastq_linux --nthreads 16 HG00514_phased_possorted_bam.bam /DATA/users/m.magnitov/hap_phen/10X/fastq/HG00514

# Puerto Rican Trio LCLs from 10X (HG00731, HG00732, HG00733):
wget -O HG00731_phased_possorted_bam.bam https://sra-pub-src-1.s3.amazonaws.com/ERR1868753/HG00731_phased_possorted_bam.bam.1
wget -O HG00732_phased_possorted_bam.bam https://sra-pub-src-1.s3.amazonaws.com/ERR1868754/HG00732_phased_possorted_bam.bam.1
wget -O HG00733_phased_possorted_bam.bam https://sra-pub-src-1.s3.amazonaws.com/ERR1868755/HG00733_phased_possorted_bam.bam.1
/DATA/users/m.magnitov/software/bamtofastq_linux --nthreads 16 HG00731_phased_possorted_bam.bam /DATA/users/m.magnitov/hap_phen/10X/fastq/HG00731
/DATA/users/m.magnitov/software/bamtofastq_linux --nthreads 16 HG00732_phased_possorted_bam.bam /DATA/users/m.magnitov/hap_phen/10X/fastq/HG00732
/DATA/users/m.magnitov/software/bamtofastq_linux --nthreads 16 HG00733_phased_possorted_bam.bam /DATA/users/m.magnitov/hap_phen/10X/fastq/HG00733

# NA20509 LCL
wget -O NA20509_phased_possorted_bam.bam http://ftp.sra.ebi.ac.uk/vol1/run/ERR259/ERR2595683/HM3CJBBXX_NA20509_longranger.bam
/DATA/users/m.magnitov/software/bamtofastq_linux --nthreads 16 NA20509_phased_possorted_bam.bam /DATA/users/m.magnitov/hap_phen/10X/fastq/NA20509

# HG03486
mkdir -p /DATA/users/m.magnitov/hap_phen/10X/fastq/HG03486
cd /DATA/users/m.magnitov/hap_phen/10X/fastq/HG03486
wget -O HG03486_S1_L001_R1_001.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/054/SRR10392754/SRR10392754_1.fastq.gz
wget -O HG03486_S1_L001_R2_001.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/054/SRR10392754/SRR10392754_2.fastq.gz
wget -O HG03486_S1_L002_R1_001.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/053/SRR10392753/SRR10392753_1.fastq.gz
wget -O HG03486_S1_L002_R2_001.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/053/SRR10392753/SRR10392753_2.fastq.gz
