# Based on the script from Wang et al., Genome Research (2019)

INPUT_PATH=/DATA/users/m.magnitov/hap_phen/simulations/linked_reads/fastq
OUTPUT_PATH=/DATA/users/m.magnitov/hap_phen/simulations/linked_reads/stLFR
GENOME=/DATA/users/m.magnitov/genomes/refdata-hg38-2.1.0/fasta/genome.fa
PICARD=/DATA/users/m.magnitov/software/picard_2.27.4.jar

mkdir -p ${OUTPUT_PATH}

#### mapping with bwa
bwa mem -t 24 -R "@RG\tID:L0\tSM:L0\tPL:COMPLETE" ${GENOME} \
	${INPUT_PATH}/stLFR_NA12878_R1.fastq.gz ${INPUT_PATH}/stLFR_NA12878_R2.fastq.gz > ${OUTPUT_PATH}/stLFR_NA12878.aln_mem.sam

perl -e '$readlen=100; $a="";$b="";$id=""; $n=0; while(<>){if(!(/^@/)){chomp;@t=split; if($t[0] ne $id){if($n>=2){print "$a\n$b\n";} $a="$_"; $n=1;}elsif(length($t[9]) == $readlen){$b="$_"; $n++;} $id=$t[0]; }else{print "$_";}} if($n>=2){print "$a\n$b\n";}' ${OUTPUT_PATH}/stLFR_NA12878.aln_mem.sam > ${OUTPUT_PATH}/stLFR_NA12878.aln.sam

#### sorting by coordinates
samtools view -buhS -t ${GENOME}.fai -T ${GENOME} ${OUTPUT_PATH}/stLFR_NA12878.aln.sam | samtools sort -@ 16 -m 5000000000 -T stLFR.sort -o ${OUTPUT_PATH}/stLFR_NA12878.sorted.bam -

#### removing duplicates in sorted bam file with picard 
java -jar ${PICARD} MarkDuplicates \
	I=${OUTPUT_PATH}/stLFR_NA12878.sorted.bam \
	O=${OUTPUT_PATH}/stLFR_NA12878.sorted.rmdup.bam \
	M=${OUTPUT_PATH}/dedup_metrics.txt \
	REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT TMP_DIR=./ \
	READ_NAME_REGEX='[a-zA-Z0-9]+_([0-9]+)#([0-9]+)_([0-9]+)_([0-9]+)'

#### filtering reads with unassigned barcodes
samtools view -@ 24 -h ${OUTPUT_PATH}/stLFR_NA12878.sorted.rmdup.bam | \
	awk -F $'\t' '($1!~/#0_0_0$/) { print }' > ${OUTPUT_PATH}/stLFR_NA12878.sorted.rmdup.rmbx.sam

#### modify bam file to add BX field with barcode information (to have a similar format as 10X data)
awk -F $'\t' '{ if ($1~/#/) {barcode=split($1,a,"#"); print $0,"BX:Z:"a[2]} else {print} }' OFS="\t" ${OUTPUT_PATH}/stLFR_NA12878.sorted.rmdup.rmbx.sam | \
	samtools view -@ 24 -h -S -b - > ${OUTPUT_PATH}/stLFR_NA12878.sorted.rmdup.rmbx.addbx.bam
samtools index ${OUTPUT_PATH}/stLFR_NA12878.sorted.rmdup.rmbx.addbx.bam

#### remove intermediate files
rm ${OUTPUT_PATH}/stLFR_NA12878.aln_mem.sam ${OUTPUT_PATH}/stLFR_NA12878.aln.sam ${OUTPUT_PATH}/stLFR_NA12878.sorted.bam ${OUTPUT_PATH}/stLFR_NA12878.sorted.rmdup.bam ${OUTPUT_PATH}/stLFR_NA12878.sorted.rmdup.rmbx.sam 
