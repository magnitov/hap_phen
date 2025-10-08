library(DESeq2)
library(ggplot2)

# Read counts and sample information
counts_allelic <- as.matrix(read.csv("/DATA/users/m.magnitov/hap_phen/TTseq/counts/counts_allelic_all_samples.txt", sep = "\t", row.names = 1))
samples_allelic <- read.csv("/DATA/users/m.magnitov/hap_phen/TTseq/ase/samples_all_allelic.txt", sep = "\t", row.names = 1)
samples_allelic$replicate <- factor(samples_allelic$replicate)
samples_allelic$allele <- factor(samples_allelic$allele)

# Create a matrix of allelic counts
dds_allelic <- DESeqDataSetFromMatrix(countData = counts_allelic,
                                      colData = samples_allelic,
                                      design = ~ allele)

# PCA of allelic counts
vsd <- vst(dds_allelic)
plotPCA(vsd, intgroup=c("replicate", "allele"), ntop = 100) +
  geom_text(aes(label=gsub(sample, "", gsub("_", "\n", name))), position = position_nudge(y = 0.1)) +
  ggtitle("PCA of haplotype counts for TT-seq samples")

# DEGs between cell lines
counts_total <- as.matrix(read.csv(sprintf("/DATA/users/m.magnitov/hap_phen/TTseq/counts/counts_total_all_samples.txt"), sep = "\t", row.names = 1))
samples_total <- read.csv(sprintf("/DATA/users/m.magnitov/hap_phen/TTseq/ase/samples_all_total.txt"), sep = "\t", row.names = 1)
samples_total$individual <- factor(samples_total$individual)
samples_total$replicate <- factor(samples_total$replicate)

dds_total <- DESeqDataSetFromMatrix(countData = counts_total,
                                    colData = samples_total,
                                    design = ~ individual)
keep <- rowSums(counts(dds_total) >= 10) >= 10
dds_total <- dds_total[keep, ]

dds_total <- estimateSizeFactors(dds_total)
dds_total <- estimateDispersions(dds_total)
dds_total <- nbinomWaldTest(dds_total)

res_NA18983_NA12878 <- results(dds_total, contrast=c("individual","NA18983","NA12878"), alpha = 0.05)
res_HG01241_NA12878 <- results(dds_total, contrast=c("individual","HG01241","NA12878"), alpha = 0.05)
res_HG02601_NA12878 <- results(dds_total, contrast=c("individual","HG02601","NA12878"), alpha = 0.05)
res_HG03464_NA12878 <- results(dds_total, contrast=c("individual","HG03464","NA12878"), alpha = 0.05)

res_HG01241_NA18983 <- results(dds_total, contrast=c("individual","HG01241","NA18983"), alpha = 0.05)
res_HG02601_NA18983 <- results(dds_total, contrast=c("individual","HG02601","NA18983"), alpha = 0.05)
res_HG03464_NA18983 <- results(dds_total, contrast=c("individual","HG03464","NA18983"), alpha = 0.05)

res_HG02601_HG01241 <- results(dds_total, contrast=c("individual","HG02601","HG01241"), alpha = 0.05)
res_HG03464_HG01241 <- results(dds_total, contrast=c("individual","HG03464","HG01241"), alpha = 0.05)

res_HG03464_HG02601 <- results(dds_total, contrast=c("individual","HG03464","HG02601"), alpha = 0.05)

write.table(as.data.frame(res_NA18983_NA12878), file = "/DATA/users/m.magnitov/hap_phen/TTseq/ase/degs_NA18983_NA12878.csv")
write.table(as.data.frame(res_HG01241_NA12878), file = "/DATA/users/m.magnitov/hap_phen/TTseq/ase/degs_HG01241_NA12878.csv")
write.table(as.data.frame(res_HG02601_NA12878), file = "/DATA/users/m.magnitov/hap_phen/TTseq/ase/degs_HG02601_NA12878.csv")
write.table(as.data.frame(res_HG03464_NA12878), file = "/DATA/users/m.magnitov/hap_phen/TTseq/ase/degs_HG03464_NA12878.csv")

write.table(as.data.frame(res_HG01241_NA18983), file = "/DATA/users/m.magnitov/hap_phen/TTseq/ase/degs_HG01241_NA18983.csv")
write.table(as.data.frame(res_HG02601_NA18983), file = "/DATA/users/m.magnitov/hap_phen/TTseq/ase/degs_HG02601_NA18983.csv")
write.table(as.data.frame(res_HG03464_NA18983), file = "/DATA/users/m.magnitov/hap_phen/TTseq/ase/degs_HG03464_NA18983.csv")

write.table(as.data.frame(res_HG02601_HG01241), file = "/DATA/users/m.magnitov/hap_phen/TTseq/ase/degs_HG02601_HG01241.csv")
write.table(as.data.frame(res_HG03464_HG01241), file = "/DATA/users/m.magnitov/hap_phen/TTseq/ase/degs_HG03464_HG01241.csv")

write.table(as.data.frame(res_HG03464_HG02601), file = "/DATA/users/m.magnitov/hap_phen/TTseq/ase/degs_HG03464_HG02601.csv")