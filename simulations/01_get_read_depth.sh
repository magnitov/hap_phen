DATA_PATH=/DATA/users/m.magnitov/hap_phen/simulations
SAMTOOLS=/DATA/users/m.magnitov/software/samtools-1.15/samtools

mkdir -p ${DATA_PATH}/depth

${SAMTOOLS} flagstat -@ 64 ${DATA_PATH}/NA12878_HiC.bam > ${DATA_PATH}/depth/flagstat_NA12878_HiC.txt
${SAMTOOLS} flagstat -@ 64 ${DATA_PATH}/NA12878_10X.bam > ${DATA_PATH}/depth/flagstat_NA12878_10X.txt
${SAMTOOLS} flagstat -@ 64 ${DATA_PATH}/NA12878_TELLseq.bam > ${DATA_PATH}/depth/flagstat_NA12878_TELLseq.txt
${SAMTOOLS} flagstat -@ 64 ${DATA_PATH}/NA12878_stLFR.bam > ${DATA_PATH}/depth/flagstat_NA12878_stLFR.txt

${SAMTOOLS} coverage --ff UNMAP,SECONDARY,SUPPLEMENTARY ${DATA_PATH}/NA12878_HiC.bam > ${DATA_PATH}/depth/coverage_NA12878_HiC.txt
${SAMTOOLS} coverage --ff UNMAP,SECONDARY,SUPPLEMENTARY ${DATA_PATH}/NA12878_10X.bam > ${DATA_PATH}/depth/coverage_NA12878_10X.txt
${SAMTOOLS} coverage --ff UNMAP,SECONDARY,SUPPLEMENTARY ${DATA_PATH}/NA12878_TELLseq.bam > ${DATA_PATH}/depth/coverage_NA12878_TELLseq.txt
${SAMTOOLS} coverage --ff UNMAP,SECONDARY,SUPPLEMENTARY ${DATA_PATH}/NA12878_stLFR.bam > ${DATA_PATH}/depth/coverage_NA12878_stLFR.txt
