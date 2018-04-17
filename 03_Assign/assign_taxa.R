
path <-  file.path("..", "00_Data") #Path to root of all data
TnFpath <- file.path(path, "clean_and_matched_1") #Path to directory of sequences that have been trimmed and filtered

library(dada2); packageVersion("dada2")

nochim <- readRDS(file.path(TnFpath, "nochim.rds"))

taxa <- assignTaxonomy(nochim, file.path(path, "silva_nr_v132_train_set.fa.gz", multithread=TRUE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
print(taxa.print)

