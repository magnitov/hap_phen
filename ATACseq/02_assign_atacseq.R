library(here)
library(data.table)
library(stringr)
library(rtracklayer)
library(GenomicAlignments)
setDTthreads(32)

assignAlleles <- function(b1, b2, outdir, outprefix) {
  
  if (!dir.exists(outdir)) dir.create(outdir)
  out01 <- here(outdir, paste(outprefix, "unassigned_hap1.bam", sep = "_"))
  out02 <- here(outdir, paste(outprefix, "unassigned_hap2.bam", sep = "_"))
  out1 <- here(outdir, paste(outprefix, "hap1.bam", sep = "_"))
  out2 <- here(outdir, paste(outprefix, "hap2.bam", sep = "_"))
  outlog <- here(outdir, paste(outprefix, ".assignment.log", sep = ""))
  
  message("Reading alignments")
  param <- ScanBamParam(
    flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isUnmappedQuery = FALSE, isDuplicate = FALSE),
    what = c("qname", "flag", "mapq", "rname", "pos", "mpos", "isize"), tag = "AS", mapqFilter = 10
  )
  message(b1)
  ga1 <- readGAlignmentPairs(b1, param = param, use.names = TRUE)
  message(b2)
  ga2 <- readGAlignmentPairs(b2, param = param, use.names = TRUE)
  
  message("Merging reads by name")
  as1 <- data.table(qname = names(ga1), AS11 = mcols(GenomicAlignments::first(ga1))$AS, AS12 = mcols(GenomicAlignments::second(ga1))$AS)
  as2 <- data.table(qname = names(ga2), AS21 = mcols(GenomicAlignments::first(ga2))$AS, AS22 = mcols(GenomicAlignments::second(ga2))$AS)
  t0 <- Sys.time()
  DT <- merge(as1, as2, by = "qname", all = TRUE)
  
  message("Assigning reads to alleles")
  DT[is.na(AS11), AS11 := 0]
  DT[is.na(AS12), AS12 := 0]
  DT[is.na(AS21), AS21 := 0]
  DT[is.na(AS22), AS22 := 0]
  DT$allele <- 0
  DT[(AS11+AS12) > (AS21+AS22), allele := 1]
  DT[(AS11+AS12) < (AS21+AS22), allele := 2]
  ga01 <- ga1[names(ga1) %in% DT[allele == 0, qname]]
  ga02 <- ga2[names(ga2) %in% DT[allele == 0, qname]]
  ga11 <- ga1[names(ga1) %in% DT[allele == 1, qname]]
  ga22 <- ga2[names(ga2) %in% DT[allele == 2, qname]]
  
  message("Writing log file")
  fwrite(
    data.table(
      c("bam1", "bam2", "total", "unassigned", "1", "2"),
      c(length(ga1), length(ga2),
        nrow(DT), sum(DT$allele == 0),
        sum(DT$allele == 1), sum(DT$allele == 2))
    ),
    outlog,
    col.names = FALSE, sep = "\t"
  )
  
  message("Writing bam files")
  export(ga01, out01, format = "bam")
  export(ga02, out02, format = "bam")
  export(ga11, out1, format = "bam")
  export(ga22, out2, format = "bam")
}


sample = "NA12878"
replicate = "rep1"
data_path = "/DATA/users/m.magnitov/hap_phen/ATACseq"

bamlist <- list.files(here(data_path, "bam"), pattern = ".bam", full.names = TRUE)
b1 <- grep(paste(sample, "_", replicate, "_", "hap1", sep = ""), bamlist, value = TRUE)
b2 <- grep(paste(sample, "_", replicate, "_", "hap2", sep = ""), bamlist, value = TRUE)
print(b1)
print(b2)

assignAlleles(b1, b2, here(data_path, "bam_assigned"), paste(sample, replicate, sep = "_"))
