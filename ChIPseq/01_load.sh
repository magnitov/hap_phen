# NA12878 ChIP-seq
mkdir -p /DATA/users/m.magnitov/hap_phen/ChIPseq/fastq/
mkdir -p /DATA/users/m.magnitov/hap_phen/ChIPseq/fastq/NA12878
cd /DATA/users/m.magnitov/hap_phen/ChIPseq/fastq/NA12878

###
### Kasowski et al. data
###
wget -O CTCF_Kasowski_rep1_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998175/SRR998175_1.fastq.gz
wget -O CTCF_Kasowski_rep1_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998175/SRR998175_2.fastq.gz
wget -O CTCF_Kasowski_rep2_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998176/SRR998176_1.fastq.gz
wget -O CTCF_Kasowski_rep2_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998176/SRR998176_2.fastq.gz
wget -O Input_Kasowski_rep1_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998196/SRR998196_1.fastq.gz
wget -O Input_Kasowski_rep1_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998196/SRR998196_2.fastq.gz
wget -O H3K27ac_Kasowski_rep1_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998177/SRR998177_1.fastq.gz
wget -O H3K27ac_Kasowski_rep1_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998177/SRR998177_2.fastq.gz
wget -O H3K27ac_Kasowski_rep2_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998178/SRR998178_1.fastq.gz
wget -O H3K27ac_Kasowski_rep2_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998178/SRR998178_2.fastq.gz
wget -O H3K27ac_Kasowski_rep3_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998179/SRR998179_1.fastq.gz
wget -O H3K27ac_Kasowski_rep3_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998179/SRR998179_2.fastq.gz
wget -O H3K27ac_Kasowski_rep4_lane1_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998180/SRR998180_1.fastq.gz
wget -O H3K27ac_Kasowski_rep4_lane1_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998180/SRR998180_2.fastq.gz
wget -O H3K27ac_Kasowski_rep4_lane2_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998181/SRR998181_1.fastq.gz
wget -O H3K27ac_Kasowski_rep4_lane2_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998181/SRR998181_2.fastq.gz
wget -O H3K27ac_Kasowski_rep4_lane3_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998182/SRR998182_1.fastq.gz
wget -O H3K27ac_Kasowski_rep4_lane3_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998182/SRR998182_2.fastq.gz
wget -O H3K27ac_Kasowski_rep4_lane4_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998183/SRR998183_1.fastq.gz
wget -O H3K27ac_Kasowski_rep4_lane4_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998183/SRR998183_2.fastq.gz
zcat H3K27ac_Kasowski_rep4_lane*_R1.fastq.gz > H3K27ac_Kasowski_rep4_R1.fastq
zcat H3K27ac_Kasowski_rep4_lane*_R2.fastq.gz > H3K27ac_Kasowski_rep4_R2.fastq
gzip H3K27ac_Kasowski_rep4_R1.fastq
gzip H3K27ac_Kasowski_rep4_R2.fastq
wget -O H3K27me3_Kasowski_rep1_lane1_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998185/SRR998185_1.fastq.gz
wget -O H3K27me3_Kasowski_rep1_lane1_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998185/SRR998185_2.fastq.gz
wget -O H3K27me3_Kasowski_rep1_lane2_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998184/SRR998184_1.fastq.gz
wget -O H3K27me3_Kasowski_rep1_lane2_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998184/SRR998184_2.fastq.gz
zcat H3K27me3_Kasowski_rep1_lane*_R1.fastq.gz > H3K27me3_Kasowski_rep1_R1.fastq
zcat H3K27me3_Kasowski_rep1_lane*_R2.fastq.gz > H3K27me3_Kasowski_rep1_R2.fastq
gzip H3K27me3_Kasowski_rep1_R1.fastq
gzip H3K27me3_Kasowski_rep1_R2.fastq
wget -O H3K27me3_Kasowski_rep2_lane1_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998187/SRR998187_1.fastq.gz
wget -O H3K27me3_Kasowski_rep2_lane1_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998187/SRR998187_2.fastq.gz
wget -O H3K27me3_Kasowski_rep2_lane2_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998186/SRR998186_1.fastq.gz
wget -O H3K27me3_Kasowski_rep2_lane2_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998186/SRR998186_2.fastq.gz
zcat H3K27me3_Kasowski_rep2_lane*_R1.fastq.gz > H3K27me3_Kasowski_rep2_R1.fastq
zcat H3K27me3_Kasowski_rep2_lane*_R2.fastq.gz > H3K27me3_Kasowski_rep2_R2.fastq
gzip H3K27me3_Kasowski_rep2_R1.fastq
gzip H3K27me3_Kasowski_rep2_R2.fastq
wget -O H3K4me3_Kasowski_rep1_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998193/SRR998193_1.fastq.gz
wget -O H3K4me3_Kasowski_rep1_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998193/SRR998193_2.fastq.gz
wget -O H3K4me3_Kasowski_rep2_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998194/SRR998194_1.fastq.gz
wget -O H3K4me3_Kasowski_rep2_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998194/SRR998194_2.fastq.gz
wget -O H3K4me3_Kasowski_rep3_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998195/SRR998195_1.fastq.gz
wget -O H3K4me3_Kasowski_rep3_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998195/SRR998195_2.fastq.gz
wget -O H3K4me1_Kasowski_rep1_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998190/SRR998190_1.fastq.gz
wget -O H3K4me1_Kasowski_rep1_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998190/SRR998190_2.fastq.gz
wget -O H3K4me1_Kasowski_rep2_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998191/SRR998191_1.fastq.gz
wget -O H3K4me1_Kasowski_rep2_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998191/SRR998191_2.fastq.gz
wget -O H3K4me1_Kasowski_rep3_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998192/SRR998192_1.fastq.gz
wget -O H3K4me1_Kasowski_rep3_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR998/SRR998192/SRR998192_2.fastq.gz

###
### Single-end datasets from ENCODE
###
wget -O BATF_ENCSR000BGT_rep1.fastq.gz https://www.encodeproject.org/files/ENCFF000NSW/@@download/ENCFF000NSW.fastq.gz
wget -O BATF_ENCSR000BGT_rep2.fastq.gz https://www.encodeproject.org/files/ENCFF000NSX/@@download/ENCFF000NSX.fastq.gz
wget -O JUNB_ENCSR897MMC_rep1.fastq.gz https://www.encodeproject.org/files/ENCFF379ZPU/@@download/ENCFF379ZPU.fastq.gz
wget -O JUNB_ENCSR897MMC_rep2.fastq.gz https://www.encodeproject.org/files/ENCFF002XBK/@@download/ENCFF002XBK.fastq.gz
wget -O SPI1_ENCSR000BGQ_rep2.fastq.gz https://www.encodeproject.org/files/ENCFF000OBU/@@download/ENCFF000OBU.fastq.gz
wget -O SPI1_ENCSR000BGQ_rep3.fastq.gz https://www.encodeproject.org/files/ENCFF000OBS/@@download/ENCFF000OBS.fastq.gz
wget -O IRF4_ENCSR000BGY_rep1.fastq.gz https://www.encodeproject.org/files/ENCFF000NWT/@@download/ENCFF000NWT.fastq.gz
wget -O IRF4_ENCSR000BGY_rep2.fastq.gz https://www.encodeproject.org/files/ENCFF000NWV/@@download/ENCFF000NWV.fastq.gz
wget -O RUNX3_ENCSR000BRI_rep1.fastq.gz https://www.encodeproject.org/files/ENCFF000OCO/@@download/ENCFF000OCO.fastq.gz
wget -O RUNX3_ENCSR000BRI_rep2.fastq.gz https://www.encodeproject.org/files/ENCFF000OCP/@@download/ENCFF000OCP.fastq.gz
wget -O POU2F2_ENCSR000BGP_rep1.fastq.gz https://www.encodeproject.org/files/ENCFF000OBK/@@download/ENCFF000OBK.fastq.gz
wget -O POU2F2_ENCSR000BGP_rep2.fastq.gz https://www.encodeproject.org/files/ENCFF000OBP/@@download/ENCFF000OBP.fastq.gz
wget -O POU2F2_ENCSR000BGP_rep3.fastq.gz https://www.encodeproject.org/files/ENCFF000OBH/@@download/ENCFF000OBH.fastq.gz

###
### Paired-end datasets from ENCODE
###
wget -O ELF1_ENCSR841NDX_rep1_R1.fastq.gz https://www.encodeproject.org/files/ENCFF164VNM/@@download/ENCFF164VNM.fastq.gz
wget -O ELF1_ENCSR841NDX_rep1_R2.fastq.gz https://www.encodeproject.org/files/ENCFF825EGN/@@download/ENCFF825EGN.fastq.gz
wget -O ELF1_ENCSR841NDX_rep2_R1.fastq.gz https://www.encodeproject.org/files/ENCFF862VDD/@@download/ENCFF862VDD.fastq.gz
wget -O ELF1_ENCSR841NDX_rep2_R2.fastq.gz https://www.encodeproject.org/files/ENCFF814SDG/@@download/ENCFF814SDG.fastq.gz

###
### Input datasets from ENCODE
###
wget -O Input_ENCSR000BGH_rep1.fastq.gz https://www.encodeproject.org/files/ENCFF000OCW/@@download/ENCFF000OCW.fastq.gz
wget -O Input_ENCSR000BGH_rep2.fastq.gz https://www.encodeproject.org/files/ENCFF000OCU/@@download/ENCFF000OCU.fastq.gz
wget -O Input_ENCSR000BGH_rep3.fastq.gz https://www.encodeproject.org/files/ENCFF000OCS/@@download/ENCFF000OCS.fastq.gz
wget -O Input_ENCSR136WQZ_rep1.fastq.gz https://www.encodeproject.org/files/ENCFF598SBH/@@download/ENCFF598SBH.fastq.gz
wget -O Input_ENCSR890ZMI_rep1.fastq.gz https://www.encodeproject.org/files/ENCFF180PND/@@download/ENCFF180PND.fastq.gz
wget -O Input_ENCSR000BVP_rep1.fastq.gz https://www.encodeproject.org/files/ENCFF000ODV/@@download/ENCFF000ODV.fastq.gz
wget -O Input_ENCSR000BVP_rep3.fastq.gz https://www.encodeproject.org/files/ENCFF000ODZ/@@download/ENCFF000ODZ.fastq.gz
wget -O Input_ENCSR000BVP_rep4.fastq.gz https://www.encodeproject.org/files/ENCFF000OEB/@@download/ENCFF000OEB.fastq.gz
wget -O Input_ENCSR956WYO_rep1_R1.fastq.gz https://www.encodeproject.org/files/ENCFF202OMW/@@download/ENCFF202OMW.fastq.gz
wget -O Input_ENCSR956WYO_rep1_R2.fastq.gz https://www.encodeproject.org/files/ENCFF349FFX/@@download/ENCFF349FFX.fastq.gz
wget -O Input_ENCSR956WYO_rep2_R1.fastq.gz https://www.encodeproject.org/files/ENCFF829GAK/@@download/ENCFF829GAK.fastq.gz
wget -O Input_ENCSR956WYO_rep2_R2.fastq.gz https://www.encodeproject.org/files/ENCFF647DBV/@@download/ENCFF647DBV.fastq.gz

### Other datasets
wget -O IRF4_MM1S_Loven_2013.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR942/SRR942953/SRR942953.fastq.gz
wget -O RELA_Zhao_rep1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/009/SRR5579179/SRR5579179.fastq.gz
wget -O RELB_Zhao_rep1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/000/SRR5579180/SRR5579180.fastq.gz
wget -O RELB_Zhao_rep2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/001/SRR5579181/SRR5579181.fastq.gz
