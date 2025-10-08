# NA12878 from Rao et al. (DpnII):
mkdir -p /DATA/users/m.magnitov/hap_phen/HiC/fastq/NA12878
cd /DATA/users/m.magnitov/hap_phen/HiC/fastq/NA12878
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1658644/SRR1658644
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1658645/SRR1658645
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1658646/SRR1658646
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1658647/SRR1658647
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1658648/SRR1658648
fastq-dump --split-files -F --gzip ./SRR1658644
fastq-dump --split-files -F --gzip ./SRR1658645
fastq-dump --split-files -F --gzip ./SRR1658646
fastq-dump --split-files -F --gzip ./SRR1658647
fastq-dump --split-files -F --gzip ./SRR1658648
zcat SRR*_1.fastq.gz > NA12878_R1.fastq
zcat SRR*_2.fastq.gz > NA12878_R2.fastq
gzip NA12878_R*
rm SRR*

# NA19238 LCLs from Gorkin et. al (HindIII)
mkdir -p /DATA/users/m.magnitov/hap_phen/HiC/fastq/NA19238
cd /DATA/users/m.magnitov/hap_phen/HiC/fastq/NA19238
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/029/SRR11935529/SRR11935529_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/029/SRR11935529/SRR11935529_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/031/SRR11935531/SRR11935531_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/031/SRR11935531/SRR11935531_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/032/SRR11935532/SRR11935532_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/032/SRR11935532/SRR11935532_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/034/SRR11935534/SRR11935534_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/034/SRR11935534/SRR11935534_2.fastq.gz
zcat SRR*_1.fastq.gz > NA19238_R1.fastq
zcat SRR*_2.fastq.gz > NA19238_R2.fastq
gzip NA19238_R*
rm SRR*

# NA19239 LCLs from Gorkin et. al (HindIII)
mkdir -p /DATA/users/m.magnitov/hap_phen/HiC/fastq/NA19239
cd /DATA/users/m.magnitov/hap_phen/HiC/fastq/NA19239
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/035/SRR11935535/SRR11935535_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/035/SRR11935535/SRR11935535_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/036/SRR11935536/SRR11935536_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/036/SRR11935536/SRR11935536_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/037/SRR11935537/SRR11935537_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/037/SRR11935537/SRR11935537_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/038/SRR11935538/SRR11935538_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/038/SRR11935538/SRR11935538_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/040/SRR11935540/SRR11935540_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/040/SRR11935540/SRR11935540_2.fastq.gz
zcat SRR*_1.fastq.gz > NA19239_R1.fastq
zcat SRR*_2.fastq.gz > NA19239_R2.fastq
gzip NA19239_R*
rm SRR*

## NA19240 LCLs from Gorkin et. al (HindIII)
mkdir -p /DATA/users/m.magnitov/hap_phen/HiC/fastq/NA19240
cd /DATA/users/m.magnitov/hap_phen/HiC/fastq/NA19240
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/042/SRR11935542/SRR11935542_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/042/SRR11935542/SRR11935542_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/043/SRR11935543/SRR11935543_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/043/SRR11935543/SRR11935543_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/044/SRR11935544/SRR11935544_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/044/SRR11935544/SRR11935544_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/046/SRR11935546/SRR11935546_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/046/SRR11935546/SRR11935546_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/047/SRR11935547/SRR11935547_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/047/SRR11935547/SRR11935547_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/048/SRR11935548/SRR11935548_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/048/SRR11935548/SRR11935548_2.fastq.gz
zcat SRR*_1.fastq.gz > NA19240_R1.fastq
zcat SRR*_2.fastq.gz > NA19240_R2.fastq
gzip NA19240_R*
rm SRR*

## HG00512 LCLs from Gorkin et. al (HindIII)
mkdir -p /DATA/users/m.magnitov/hap_phen/HiC/fastq/HG00512
cd /DATA/users/m.magnitov/hap_phen/HiC/fastq/HG00512
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/009/SRR8768099/SRR8768099_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/009/SRR8768099/SRR8768099_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/000/SRR8768100/SRR8768100_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/000/SRR8768100/SRR8768100_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/001/SRR8768101/SRR8768101_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/001/SRR8768101/SRR8768101_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/002/SRR8768102/SRR8768102_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/002/SRR8768102/SRR8768102_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/003/SRR8768103/SRR8768103_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/003/SRR8768103/SRR8768103_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/004/SRR8768104/SRR8768104_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/004/SRR8768104/SRR8768104_2.fastq.gz
zcat SRR*_1.fastq.gz > HG00512_R1.fastq
zcat SRR*_2.fastq.gz > HG00512_R2.fastq
gzip HG00512_R*
rm SRR*

## HG00513 LCLs from Gorkin et. al (HindIII)
mkdir -p /DATA/users/m.magnitov/hap_phen/HiC/fastq/HG00513
cd /DATA/users/m.magnitov/hap_phen/HiC/fastq/HG00513
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/005/SRR8768105/SRR8768105_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/005/SRR8768105/SRR8768105_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/006/SRR8768106/SRR8768106_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/006/SRR8768106/SRR8768106_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/007/SRR8768107/SRR8768107_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/007/SRR8768107/SRR8768107_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/008/SRR8768108/SRR8768108_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/008/SRR8768108/SRR8768108_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/009/SRR8768109/SRR8768109_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/009/SRR8768109/SRR8768109_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/000/SRR8768110/SRR8768110_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/000/SRR8768110/SRR8768110_2.fastq.gz
zcat SRR*_1.fastq.gz > HG00513_R1.fastq
zcat SRR*_2.fastq.gz > HG00513_R2.fastq
gzip HG00513_R*
rm SRR*

## HG00514 LCLs from Gorkin et. al (HindIII)
mkdir -p /DATA/users/m.magnitov/hap_phen/HiC/fastq/HG00514
cd /DATA/users/m.magnitov/hap_phen/HiC/fastq/HG00514
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/001/SRR8768111/SRR8768111_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/001/SRR8768111/SRR8768111_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/002/SRR8768112/SRR8768112_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/002/SRR8768112/SRR8768112_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/003/SRR8768113/SRR8768113_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/003/SRR8768113/SRR8768113_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/004/SRR8768114/SRR8768114_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/004/SRR8768114/SRR8768114_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/005/SRR8768115/SRR8768115_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/005/SRR8768115/SRR8768115_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/006/SRR8768116/SRR8768116_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/006/SRR8768116/SRR8768116_2.fastq.gz
zcat SRR*_1.fastq.gz > HG00514_R1.fastq
zcat SRR*_2.fastq.gz > HG00514_R2.fastq
gzip HG00514_R*
rm SRR*

## HG00731 LCLs from Gorkin et. al (HindIII)
mkdir -p /DATA/users/m.magnitov/hap_phen/HiC/fastq/HG00731
cd /DATA/users/m.magnitov/hap_phen/HiC/fastq/HG00731
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/007/SRR8768117/SRR8768117_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/007/SRR8768117/SRR8768117_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/008/SRR8768118/SRR8768118_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/008/SRR8768118/SRR8768118_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/009/SRR8768119/SRR8768119_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/009/SRR8768119/SRR8768119_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/000/SRR8768120/SRR8768120_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/000/SRR8768120/SRR8768120_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/001/SRR8768121/SRR8768121_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/001/SRR8768121/SRR8768121_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/002/SRR8768122/SRR8768122_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/002/SRR8768122/SRR8768122_2.fastq.gz
zcat SRR*_1.fastq.gz > HG00731_R1.fastq
zcat SRR*_2.fastq.gz > HG00731_R2.fastq
gzip HG00731_R*
rm SRR*

## HG00732 LCLs from Gorkin et. al (HindIII)
mkdir -p /DATA/users/m.magnitov/hap_phen/HiC/fastq/HG00732
cd /DATA/users/m.magnitov/hap_phen/HiC/fastq/HG00732
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/003/SRR8768123/SRR8768123_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/003/SRR8768123/SRR8768123_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/004/SRR8768124/SRR8768124_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/004/SRR8768124/SRR8768124_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/005/SRR8768125/SRR8768125_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/005/SRR8768125/SRR8768125_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/006/SRR8768126/SRR8768126_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/006/SRR8768126/SRR8768126_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/007/SRR8768127/SRR8768127_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/007/SRR8768127/SRR8768127_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/008/SRR8768128/SRR8768128_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/008/SRR8768128/SRR8768128_2.fastq.gz
zcat SRR*_1.fastq.gz > HG00732_R1.fastq
zcat SRR*_2.fastq.gz > HG00732_R2.fastq
gzip HG00732_R*
rm SRR*

## HG00733 LCLs from Gorkin et. al (HindIII)
mkdir -p /DATA/users/m.magnitov/hap_phen/HiC/fastq/HG00733
cd /DATA/users/m.magnitov/hap_phen/HiC/fastq/HG00733
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/009/SRR8768129/SRR8768129_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/009/SRR8768129/SRR8768129_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/000/SRR8768130/SRR8768130_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/000/SRR8768130/SRR8768130_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/001/SRR8768131/SRR8768131_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/001/SRR8768131/SRR8768131_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/002/SRR8768132/SRR8768132_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/002/SRR8768132/SRR8768132_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/003/SRR8768133/SRR8768133_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/003/SRR8768133/SRR8768133_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/004/SRR8768134/SRR8768134_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR876/004/SRR8768134/SRR8768134_2.fastq.gz
zcat SRR*_1.fastq.gz > HG00733_R1.fastq
zcat SRR*_2.fastq.gz > HG00733_R2.fastq
gzip HG00733_R*
rm SRR*

# NA20509 
mkdir -p /DATA/users/m.magnitov/hap_phen/HiC/fastq/NA20509
cd /DATA/users/m.magnitov/hap_phen/HiC/fastq/NA20509
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/ERR4501132/ERR4501132
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/ERR4471592/ERR4471592
fastq-dump --split-files -F --gzip ./ERR4501132
fastq-dump --split-files -F --gzip ./ERR4471592
mv ERR4501132_1.fastq.gz NA20509_R1.fastq.gz
mv ERR4471592_1.fastq.gz NA20509_R2.fastq.gz

# HG03486
mkdir -p /DATA/users/m.magnitov/hap_phen/HiC/fastq/HG03486
cd /DATA/users/m.magnitov/hap_phen/HiC/fastq/HG03486
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR117/029/SRR11777629/SRR11777629_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR117/029/SRR11777629/SRR11777629_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR117/028/SRR11777628/SRR11777628_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR117/028/SRR11777628/SRR11777628_2.fastq.gz
zcat SRR*_1.fastq.gz > HG03486_R1.fastq
zcat SRR*_2.fastq.gz > HG03486_R2.fastq
gzip HG03486_R*
rm SRR*

