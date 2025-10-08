#!/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/personal_genomes
GENOME=/DATA/users/m.magnitov/genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
MINIMAP2=/DATA/users/m.magnitov/software/minimap2-2.24_x64-linux/minimap2
PAF2CHAIN=/home/m.magnitov/.cargo/bin/paf2chain

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
    # Align haplotypes to reference genome with minimap2
    ${MINIMAP2} -t 32 -cx asm5 --cs ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap1.fa ${GENOME} > ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap1.paf
    ${MINIMAP2} -t 32 -cx asm5 --cs ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap2.fa ${GENOME} > ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap2.paf
    
    # Convert obtained PAF files to CHAIN files compatible with liftOver
    ${PAF2CHAIN} -i ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap1.paf > ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap1_to_ref.chain
    ${PAF2CHAIN} -i ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap2.paf > ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap2_to_ref.chain
    awk '{ if ($0 ~ /chain/) print $1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "$10" "$11" "$12" "$13; else print $0 }' \
	    ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap1_to_ref.chain > ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap1_to_ref.modified.chain
    awk '{ if ($0 ~ /chain/) print $1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "$10" "$11" "$12" "$13; else print $0 }' \
	    ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap2_to_ref.chain > ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap2_to_ref.modified.chain
    mv ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap1_to_ref.modified.chain ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap1_to_ref.chain
    mv ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap2_to_ref.modified.chain ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap2_to_ref.chain
done
