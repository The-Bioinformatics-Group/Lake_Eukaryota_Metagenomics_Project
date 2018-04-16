path <-  file.path("..", "00_Data") #Path to root of all data
TnFpath1 <- file.path(path, "clean_and_matched_1") #Path to directory of sequences that have been trimmed and filtered
TnFpath2 <- file.path(path, "clean_and_matched_2")

library(dada2); packageVersion("dada2")

dadaF1 <- readRDS(file.path(TnFpath1, "dadaF.rds"))
mergers1 <- readRDS(file.path(TnFpath1, "mergers.rds"))
seqtab1 <- readRDS(file.path(TnFpath1, "seqTab.rds"))

dadaF2 <- readRDS(file.path(TnFpath2, "dadaF.rds"))
mergers2 <- readRDS(file.path(TnFpath2, "mergers.rds"))
seqtab2 <- readRDS(file.path(TnFpath2, "seqTab.rds"))

print(Analysis alt 1)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaF1, getN), sapply(mergers1, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

print(Analysis alt 2)
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)