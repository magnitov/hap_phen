library(here)
library(data.table)
library(stringr)
library(rtracklayer)
library(GenomicAlignments)
setDTthreads(32)

assignAlleles <- function(b1, b2, outdir, outprefix) {
  
  if (!dir.exists(outdir)) dir.create(outdir)
  #out01 <- here(outdir, paste(outprefix, "unassigned_hap1.bam", sep = "_"))
  #out02 <- here(outdir, paste(outprefix, "unassigned_hap2.bam", sep = "_"))
  out1 <- here(outdir, paste(outprefix, "hap1.bam", sep = "_"))
  out2 <- here(outdir, paste(outprefix, "hap2.bam", sep = "_"))
  outlog <- here(outdir, paste(outprefix, ".assignment.log", sep = ""))
  
  message("Reading alignments")
  param <- ScanBamParam(
    flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE),
    what = c("qname", "flag", "mapq", "rname", "pos", "mpos", "isize"), tag = "AS", mapqFilter = 10
  )
  message(b1)
  ga1 <- readGAlignments(b1, param = param, use.names = TRUE)
  message(b2)
  ga2 <- readGAlignments(b2, param = param, use.names = TRUE)
  
  message("Merging reads by name")
  as1 <- data.table(qname = names(ga1), AS1 = mcols(ga1)$AS)
  as2 <- data.table(qname = names(ga2), AS2 = mcols(ga2)$AS)
  t0 <- Sys.time()
  DT <- merge(as1, as2, by = "qname", all = TRUE)
  print(length(DT))
  
  message("Assigning reads to alleles")
  DT[is.na(AS1), AS1 := 0]
  DT[is.na(AS2), AS2 := 0]
  DT$allele <- 0
  DT[AS1 > AS2, allele := 1]
  DT[AS1 < AS2, allele := 2]
  #ga01 <- ga1[names(ga1) %in% DT[allele == 0, qname]]
  #ga02 <- ga2[names(ga2) %in% DT[allele == 0, qname]]
  ga11 <- ga1[names(ga1) %in% DT[allele == 1, qname]]
  ga22 <- ga2[names(ga2) %in% DT[allele == 2, qname]]
  print(length(ga11))
  print(length(ga22))
  
  message("Writing bam files")
  #export(ga01, out01, format = "bam")
  #export(ga02, out02, format = "bam")
  export(ga11, out1, format = "bam")
  export(ga22, out2, format = "bam")
  
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
}


data_path = "/DATA/users/m.magnitov/hap_phen/ChIPseq"

for (sample in c("RELA_Zhao_rep1", "RELB_Zhao_rep1", "RELB_Zhao_rep2", "IRF4_MM1S_Loven_2013", "BATF_ENCSR000BGT_rep1", "BATF_ENCSR000BGT_rep2", "IRF4_ENCSR000BGY_rep1", "IRF4_ENCSR000BGY_rep2", "JUNB_ENCSR897MMC_rep1", "JUNB_ENCSR897MMC_rep2", "POU2F2_ENCSR000BGP_rep1", "POU2F2_ENCSR000BGP_rep2", "POU2F2_ENCSR000BGP_rep3", "RUNX3_ENCSR000BRI_rep1", "RUNX3_ENCSR000BRI_rep2", "SPI1_ENCSR000BGQ_rep2", "SPI1_ENCSR000BGQ_rep3", "Input_ENCSR000BGH_rep1", "Input_ENCSR000BGH_rep2", "Input_ENCSR000BGH_rep3", "Input_ENCSR000BVP_rep1", "Input_ENCSR000BVP_rep3", "Input_ENCSR000BVP_rep4", "Input_ENCSR136WQZ_rep1", "Input_ENCSR890ZMI_rep1")) {
    bamlist <- list.files(here(data_path, "bam"), pattern = ".bam", full.names = TRUE)
    b1 <- grep(paste(sample, "_", "hap1", sep = ""), bamlist, value = TRUE)
    b2 <- grep(paste(sample, "_", "hap2", sep = ""), bamlist, value = TRUE)
    print(b1)
    print(b2)

    assignAlleles(b1, b2, here(data_path, "bam_assigned"), paste(sample, sep = "_"))
}
