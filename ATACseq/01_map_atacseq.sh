SAMPLE=$1
REPLICATE=$2
INDIVIDUAL=$3
DATA_PATH=/DATA/users/m.magnitov/hap_phen/ATACseq
GENOME_PATH=/DATA/users/m.magnitov/hap_phen/personal_genomes
REFERENCE_GENOME=/DATA/users/m.magnitov/genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
PICARD=/DATA/users/m.magnitov/software/picard_2.27.4.jar

mkdir -p ${DATA_PATH}/bam
mkdir -p ${DATA_PATH}/bam_assigned

# Mapping reads
bwa mem -t 16 -M ${REFERENCE_GENOME} \
	${DATA_PATH}/fastq/${INDIVIDUAL}/${SAMPLE}_${REPLICATE}_R1.fastq.gz ${DATA_PATH}/fastq/${INDIVIDUAL}/${SAMPLE}_${REPLICATE}_R2.fastq.gz |\
	awk '{ if ($3!="chrM") print $0}' |\
	samtools view -@ 16 -h -S -b -f 2 -q 10 |\
	samtools sort -@ 16 - > ${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_ref_Aligned.sortedByCoord.mapq.bam

bwa mem -t 16 -M ${GENOME_PATH}/${INDIVIDUAL}/genome/${INDIVIDUAL}_hap1.fa \
	${DATA_PATH}/fastq/${INDIVIDUAL}/${SAMPLE}_${REPLICATE}_R1.fastq.gz ${DATA_PATH}/fastq/${INDIVIDUAL}/${SAMPLE}_${REPLICATE}_R2.fastq.gz |\
	awk '{ if ($3!="chrM") print $0}' |\
	samtools view -@ 16 -h -S -b -f 2 -q 10 |\
	samtools sort -@ 16 - > ${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_hap1_Aligned.sortedByCoord.mapq.bam

bwa mem -t 16 -M ${GENOME_PATH}/${INDIVIDUAL}/genome/${INDIVIDUAL}_hap2.fa \
	${DATA_PATH}/fastq/${INDIVIDUAL}/${SAMPLE}_${REPLICATE}_R1.fastq.gz ${DATA_PATH}/fastq/${INDIVIDUAL}/${SAMPLE}_${REPLICATE}_R2.fastq.gz |\
	awk '{ if ($3!="chrM") print $0}' |\
	samtools view -@ 16 -h -S -b -f 2 -q 10 |\
	samtools sort -@ 16 - > ${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_hap2_Aligned.sortedByCoord.mapq.bam
    
# Mark duplicated read pairs
java -jar ${PICARD} MarkDuplicates REMOVE_DUPLICATES=true \
	I=${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_ref_Aligned.sortedByCoord.mapq.bam \
	O=${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_ref_Aligned.sortedByCoord.mapq.nodups.bam \
	M=${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_ref_dedupMetrics.txt

java -jar ${PICARD} MarkDuplicates REMOVE_DUPLICATES=true \
	I=${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_hap1_Aligned.sortedByCoord.mapq.bam \
	O=${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_hap1_Aligned.sortedByCoord.mapq.nodups.bam \
	M=${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_hap1_dedupMetrics.txt

java -jar ${PICARD} MarkDuplicates REMOVE_DUPLICATES=true \
	I=${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_hap2_Aligned.sortedByCoord.mapq.bam \
	O=${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_hap2_Aligned.sortedByCoord.mapq.nodups.bam \
	M=${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_hap2_dedupMetrics.txt

# Remove intermediate files
rm ${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_*_Aligned.sortedByCoord.mapq.bam
