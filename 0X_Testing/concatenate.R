path <-  file.path("..", "00_Data") #Path to root of all data
TnFpath <- file.path(path, "clean_and_matched_1") #Path to directory of sequences that have been trimmed and filtered

library(dada2); packageVersion("dada2")

dadaF <- readRDS(file.path(TnFpath, "dadaF.rds"))
dadaR <- readRDS(file.path(TnFpath, "dadaR.rds"))
derepF <- readRDS(file.path(TnFpath, "derepF.rds"))
derepR <- readRDS(file.path(TnFpath, "derepR.rds"))

mergers <- mergePairs(dadaF, derepF, dadaR, derepR, justConcatenate=TRUE, verbose=TRUE)

saveRDS(mergers, file.path(TnFpath, "mergers.rds"))