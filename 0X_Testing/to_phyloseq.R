path <-  file.path("..", "00_Data") #Path to root of all data
TnFpath <- file.path(path, "clean_and_matched_2") #Path to directory of sequences that have been trimmed and filtered

library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

nochim <- readRDS(file.path(TnFpath, "nochim.rds"))
taxa <- readRDS(file.path(TnFpath, "taxa.rds"))

samples.out <- rownames(nochim)

subject <- sapply(strsplit(sample.out, "_"), `[`, 1)

samdf <- data.frame(Subject=subject)
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(nochim, taxa_are_rows=FALSE), sample_data(samdf), tax_table(taxa))