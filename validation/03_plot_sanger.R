# Based on the script developed by Moreno Martinovic

setwd("/DATA/users/m.magnitov/hap_phen/validation")

library(here)
library(sangerseqR)
library(Biostrings)
library(data.table)
library(ggplot2)

theme_set(theme_grey())
theme_update(text = element_text(colour = "black"),
             line = element_line(colour = "black"),
             axis.line  = element_line(colour = "black"),
             axis.ticks = element_line(colour = "black"),
             axis.text  = element_text(colour = "black"),
             legend.key = element_blank(),
             legend.background = element_rect(colour = NA, fill = NA),
             panel.background = element_blank(),
             panel.grid.major = element_line(colour = "grey95"),
             panel.grid.minor = element_blank(),
             plot.background  = element_blank(),
             strip.background = element_blank(),
             strip.text = element_text(colour = "black"))

# 2514ZAC159_premix.ab1 - cDNA wild-type
# 2514ZAC158_premix.ab1 - cDNA edited pool
# 2514ZAC171_premix.ab1 - gDNA wild-type
# 2514ZAC170_premix.ab1 - gDNA edited pool

## Read Sanger files
# gDNA
infiles <- list.files(".", pattern = "ab1$", full.names = TRUE)
infiles <- infiles[3:4]
names(infiles) <- c("Pool", "WT")
scf <- lapply(infiles, readsangerseq)
# cDNA
#infiles <- list.files(".", pattern = "ab1$", full.names = TRUE)
#infiles <- infiles[1:2]
#names(infiles) <- c("Pool", "WT")
#scf <- lapply(infiles, readsangerseq)

## Get traces
traces <- lapply(scf, function(x) {
  x <- traceMatrix(x)
  colnames(x) <- c("A", "C", "G", "T")
  x
})
trace_df <- lapply(traces, reshape2::melt)
trace_df <- rbindlist(trace_df, idcol = "sample")
setnames(trace_df, "Var2", "base")
trace_df[, base := c("A", "C", "G", "T")[base]]

## Get base positions
positions <- lapply(scf, function(x) data.frame(Var1 = peakPosMatrix(x)[,1]))
positions <- rbindlist(positions, idcol = "sample")

## Get base calls
bases <- lapply(scf, primarySeq, string = TRUE)
bases <- lapply(bases, function(x) unlist(strsplit(x, "")))
positions$base <- unlist(bases)

## Plot example chromatograms
ggplot(trace_df[Var1 > 6858 & Var1 < 7135],
       aes(Var1, value, colour = base)) +
  geom_line(aes(group = base)) +
  geom_text(data = positions[Var1 > 6858 & Var1 < 7135],
            aes(Var1, -200, label = base)) +
  scale_x_continuous(expand = c(0, 0), name = NULL, breaks = NULL) +
  scale_colour_manual(values = c("green3", "dodgerblue3", "black", "brown2")) +
  theme(aspect.ratio = 1/5,
        axis.line.x = element_blank(),
        panel.grid.major = element_blank()) +
  facet_grid(sample ~ .)
ggsave(here("figures", "sanger_gDNA_edited_variant_tracks.pdf"), width = 7, height = 4)

ggplot(trace_df[Var1 > 4848 & Var1 < 5050],
       aes(Var1, value, colour = base)) +
  geom_line(aes(group = base)) +
  geom_text(data = positions[Var1 > 4848 & Var1 < 5050],
            aes(Var1, -200, label = base)) +
  scale_x_continuous(expand = c(0, 0), name = NULL, breaks = NULL) +
  scale_colour_manual(values = c("green3", "dodgerblue3", "black", "brown2")) +
  theme(aspect.ratio = 1/5,
        axis.line.x = element_blank(),
        panel.grid.major = element_blank()) +
  facet_grid(sample ~ .)
ggsave(here("figures", "sanger_gDNA_non_edited_variant_tracks.pdf"), width = 7, height = 4)

## Plot example chromatograms
#ggplot(trace_df[Var1 > 898 & Var1 < 1059],
#       aes(Var1, value, colour = base)) +
#  geom_line(aes(group = base)) +
#  geom_text(data = positions[Var1 > 898 & Var1 < 1059],
#            aes(Var1, -200, label = base)) +
#  scale_x_continuous(expand = c(0, 0), name = NULL, breaks = NULL) +
#  scale_colour_manual(values = c("green3", "dodgerblue3", "black", "brown2")) +
#  theme(aspect.ratio = 1/5,
#        axis.line.x = element_blank(),
#        panel.grid.major = element_blank()) +
#  facet_grid(sample ~ .)
#ggsave(here("figures", "sanger_cDNA_tracks.pdf"), width = 7, height = 4)