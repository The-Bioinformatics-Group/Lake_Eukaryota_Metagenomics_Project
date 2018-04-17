path <-  file.path("..", "00_Data") #Path to root of all data
TnFpath <- file.path(path, "clean_and_matched_2") #Path to directory of sequences that have been trimmed and filtered

library(dada2); packageVersion("dada2")

dadaF <- readRDS(file.path(TnFpath, "dadaF.rds"))
dadaR <- readRDS(file.path(TnFpath, "dadaR.rds"))
derepF <- readRDS(file.path(TnFpath, "derepF.rds"))
derepR <- readRDS(file.path(TnFpath, "derepR.rds"))
mergers <- readRDS(file.path(TnFpath, "mergers.rds"))

seqtab <- makeSequenceTable(mergers)

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

saveRDS(seqtab, file.path(TnFpath, "seqTab.rds"))
saveRDS(seqtab.nochim, file.path(TnFpath, "nochim.rds"))
