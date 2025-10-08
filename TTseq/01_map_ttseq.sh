SAMPLE=$1
REPLICATE=$2
DATA_PATH=/DATA/users/m.magnitov/hap_phen/TTseq
GENOME_PATH=/DATA/users/m.magnitov/hap_phen/personal_genomes
REFERENCE_GENOME=/DATA/users/m.magnitov/genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set_STARgenome
STAR=/DATA/users/m.magnitov/software/STAR-2.7.10a/bin/Linux_x86_64/STAR
PICARD=/DATA/users/m.magnitov/software/picard_2.27.4.jar

mkdir -p ${DATA_PATH}/bam
mkdir -p ${DATA_PATH}/bam_assigned

# Mapping reads
${STAR} --runThreadN 32 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat \
	--genomeDir ${REFERENCE_GENOME} \
	--outFileNamePrefix ${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_ref_ \
	--readFilesIn ${DATA_PATH}/fastq/${SAMPLE}/${SAMPLE}_${REPLICATE}_R1.fastq.gz ${DATA_PATH}/fastq/${SAMPLE}/${SAMPLE}_${REPLICATE}_R2.fastq.gz \
	--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 999 --alignSJDBoverhangMin 1 \
	--outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04
    
${STAR} --runThreadN 32 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat \
	--genomeDir ${GENOME_PATH}/${SAMPLE}/genome/${SAMPLE}_hap1_STARgenome \
	--outFileNamePrefix ${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_hap1_ \
	--readFilesIn ${DATA_PATH}/fastq/${SAMPLE}/${SAMPLE}_${REPLICATE}_R1.fastq.gz ${DATA_PATH}/fastq/${SAMPLE}/${SAMPLE}_${REPLICATE}_R2.fastq.gz \
	--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 999 --alignSJDBoverhangMin 1 \
	--outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04

${STAR} --runThreadN 32 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat \
	--genomeDir ${GENOME_PATH}/${SAMPLE}/genome/${SAMPLE}_hap2_STARgenome \
	--outFileNamePrefix ${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_hap2_ \
	--readFilesIn ${DATA_PATH}/fastq/${SAMPLE}/${SAMPLE}_${REPLICATE}_R1.fastq.gz ${DATA_PATH}/fastq/${SAMPLE}/${SAMPLE}_${REPLICATE}_R2.fastq.gz \
	--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 999 --alignSJDBoverhangMin 1 \
	--outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04

# Filtering reads by mapping quality
samtools view -h -S -b -q 10 -@ 32 -o ${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_ref_Aligned.sortedByCoord.mapq.bam \
	${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_ref_Aligned.sortedByCoord.out.bam
samtools view -h -S -b -q 10 -@ 32 -o ${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_hap1_Aligned.sortedByCoord.mapq.bam \
	${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_hap1_Aligned.sortedByCoord.out.bam
samtools view -h -S -b -q 10 -@ 32 -o ${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_hap2_Aligned.sortedByCoord.mapq.bam \
	${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_hap2_Aligned.sortedByCoord.out.bam

# Mark duplicated read pairs
java -jar ${PICARD} MarkDuplicates REMOVE_DUPLICATES= true \
	I=${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_ref_Aligned.sortedByCoord.mapq.bam \
	O=${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_ref_Aligned.sortedByCoord.mapq.nodups.bam \
	M=${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_ref_dedupMetrics.txt
    
java -jar ${PICARD} MarkDuplicates REMOVE_DUPLICATES= true \
	I=${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_hap1_Aligned.sortedByCoord.mapq.bam \
	O=${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_hap1_Aligned.sortedByCoord.mapq.nodups.bam \
	M=${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_hap1_dedupMetrics.txt

java -jar ${PICARD} MarkDuplicates REMOVE_DUPLICATES= true \
	I=${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_hap2_Aligned.sortedByCoord.mapq.bam \
	O=${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_hap2_Aligned.sortedByCoord.mapq.nodups.bam \
	M=${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_hap2_dedupMetrics.txt

# Remove intermediate files
rm ${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_*_SJ.out.tab
rm ${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_*_Log.progress.out
rm ${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_*_Aligned.sortedByCoord.out.bam
rm ${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_*_Aligned.sortedByCoord.mapq.bam
