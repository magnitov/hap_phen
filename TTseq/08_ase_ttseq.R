# Based on the script developed by Robin H. van der Weide

library(DESeq2)

# Define sample
sample <- "HG03464"

# Read counts and sample information
counts_allelic <- as.matrix(read.csv(sprintf("/DATA/users/m.magnitov/hap_phen/TTseq/counts/counts_allelic_%s.txt", sample), sep = "\t", row.names = 1))
counts_total <- as.matrix(read.csv(sprintf("/DATA/users/m.magnitov/hap_phen/TTseq/counts/counts_total_%s.txt", sample), sep = "\t", row.names = 1))

samples_allelic <- read.csv(sprintf("/DATA/users/m.magnitov/hap_phen/TTseq/ase/samples_%s_allelic.txt", sample), sep = "\t", row.names = 1)
samples_allelic$replicate <- factor(samples_allelic$replicate)
samples_allelic$allele <- factor(samples_allelic$allele)

samples_total <- read.csv(sprintf("/DATA/users/m.magnitov/hap_phen/TTseq/ase/samples_%s.txt", sample), sep = "\t", row.names = 1)
samples_total$replicate <- factor(samples_total$replicate)

# Create a matrix of allelic counts and filter low counts
dds_allelic <- DESeqDataSetFromMatrix(countData = counts_allelic,
                                      colData = samples_allelic,
                                      design = ~ replicate+allele)

keep <- rowSums(counts(dds_allelic) >= 5) >= 2
dds_allelic <- dds_allelic[keep, ]

vsd <- vst(dds_allelic, blind=FALSE)
plotPCA(vsd, intgroup=c("replicate", "allele")) +
  geom_text(aes(label=gsub(sample, "", gsub("_", "\n", name))), position = position_nudge(y = 2)) +
  ggtitle(sprintf("PCA of haplotype counts for %s TT-seq samples", sample))

# Create a matrix of total counts and estimate size factors and dispersion
dds_total <- DESeqDataSetFromMatrix(countData = counts_total,
                                    colData = samples_total, 
                                    design = ~ 1)
dds_total <- estimateSizeFactors(dds_total)
dds_total <- estimateDispersions(dds_total, fitType = 'local')
dds_total <- dds_total[rownames(dds_allelic), ]

# Apply size factors and dispersion to allelic counts
sizeFactors(dds_allelic) <- sizeFactors(dds_total)[c(1, 1, 2, 2)]
rowData(dds_allelic) <- rowData(dds_total)

# Perform differential allelic expression test
dds_allelic <- nbinomWaldTest(dds_allelic)

# Extract results and save them
res <- results(dds_allelic, alpha = 0.1)
summary(res)

hist(res$pvalue, main = sprintf("Histogram of p-values for %s TT-seq", sample))

res <- res[!is.na(res$padj), ]
write.table(as.data.frame(res), file = sprintf("/DATA/users/m.magnitov/hap_phen/TTseq/ase/ase_%s.csv", sample))
