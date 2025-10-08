DATA_PATH=/DATA/users/m.magnitov/hap_phen/simulations
SAMTOOLS=/DATA/users/m.magnitov/software/samtools-1.15/samtools 

# mean depth for chr1 from samtools coverage
DEPTH_10X=54.0805
DEPTH_stLFR=38.6251
DEPTH_TELLseq=34.4986
DEPTH_HiC=35.7642

mkdir -p ${DATA_PATH}/subsample_10X
mkdir -p ${DATA_PATH}/subsample_stLFR
mkdir -p ${DATA_PATH}/subsample_TELLseq
mkdir -p ${DATA_PATH}/subsample_HiC

for SAMPLING_NUM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
	for DEPTH in 1 2 3 4 5 6 7 8 10 12 16 20 24 32
	do
		# Subsample 10X
		FRACTION_10X=$(echo "scale=10;$DEPTH/$DEPTH_10X" | bc)
		${SAMTOOLS} view -@ 32 -h -S -b -F 4 -F 256 -F 2048 \
			--subsample ${FRACTION_10X} --subsample-seed ${SAMPLING_NUM} \
			-o ${DATA_PATH}/subsample_10X/NA12878_10X_chr1_${DEPTH}x_ds${SAMPLING_NUM}.bam \
			${DATA_PATH}/NA12878_10X.bam chr1
		${SAMTOOLS} index -@ 64	${DATA_PATH}/subsample_10X/NA12878_10X_chr1_${DEPTH}x_ds${SAMPLING_NUM}.bam

		# Subsample stLFR
		FRACTION_stLFR=$(echo "scale=10;$DEPTH/$DEPTH_stLFR" | bc)
		${SAMTOOLS} view -@ 32 -h -S -b -F 4 -F 256 -F 2048 \
			--subsample ${FRACTION_stLFR} --subsample-seed ${SAMPLING_NUM} \
			-o ${DATA_PATH}/subsample_stLFR/NA12878_stLFR_chr1_${DEPTH}x_ds${SAMPLING_NUM}.bam \
			${DATA_PATH}/NA12878_stLFR.bam chr1
		${SAMTOOLS} index -@ 64 ${DATA_PATH}/subsample_stLFR/NA12878_stLFR_chr1_${DEPTH}x_ds${SAMPLING_NUM}.bam

		# Subsample TELL-seq
		FRACTION_TELLseq=$(echo "scale=10;$DEPTH/$DEPTH_TELLseq" | bc)
		${SAMTOOLS} view -@ 32 -h -S -b -F 4 -F 256 -F 2048 \
			--subsample ${FRACTION_TELLseq} --subsample-seed ${SAMPLING_NUM} \
			-o ${DATA_PATH}/subsample_TELLseq/NA12878_TELLseq_chr1_${DEPTH}x_ds${SAMPLING_NUM}.bam \
			${DATA_PATH}/NA12878_TELLseq.bam chr1
		${SAMTOOLS} index -@ 64 ${DATA_PATH}/subsample_TELLseq/NA12878_TELLseq_chr1_${DEPTH}x_ds${SAMPLING_NUM}.bam
	
		# Subsample Hi-C
		FRACTION_HiC=$(echo "scale=10;$DEPTH/$DEPTH_HiC" | bc)
		${SAMTOOLS} view -@ 32 -h -S -b -F 4 -F 256 -F 2048 \
			--subsample ${FRACTION_HiC} --subsample-seed ${SAMPLING_NUM} \
			-o ${DATA_PATH}/subsample_HiC/NA12878_HiC_chr1_${DEPTH}x_ds${SAMPLING_NUM}.bam \
			${DATA_PATH}/NA12878_HiC.bam chr1
		${SAMTOOLS} index -@ 64 ${DATA_PATH}/subsample_HiC/NA12878_HiC_chr1_${DEPTH}x_ds${SAMPLING_NUM}.bam
	done
done

