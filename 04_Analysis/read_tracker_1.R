## Prints the number of reads in each sample after each step in the data processing
## of data set 1/ trimming method 1
## Not fully implemented

path <-  file.path("..", "00_Data") #Path to root of all data
TnFpath1 <- file.path(path, "clean_and_matched_1") #Path to directory of sequences that have been trimmed and filtered

files <- list.files(TnFpath1, pattern="R1_001.fastq.gz")

library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")

#Raw

#Adapters removed
AdapFiles <- dir("../00_Data/AdaptersRemoved", "*R1_001.fastq.gz")
qaSummary <- qa(file.path(path, "AdaptersRemoved", AdapFiles), type="fastq")
NoAdap <- head(qaSummary[["readCounts"]])

#Trimmed 1
TrimFiles <- dir("../00_Data/trimmedTailsw", "*R1_001.fastq.gz")
qaSummary <- qa(file.path(path, "trimmedTailsw", TrimFiles), type="fastq")
Trim <- head(qaSummary[["readCounts"]])

#Filtered 1
Filt1Files <- dir("../00_Data/longSeq", "*R1_001.fastq.gz")
qaSummary <- qa(file.path(path, "longSeq", Filt1Files), type="fastq")
Filt1 <- head(qaSummary[["readCounts"]])

#Filtered 1
Filt2Files <- dir("../00_Data/clean_and_matched_1", "*R1_001.fastq.gz")
qaSummary <- qa(file.path(path, "clean_and_matched_1", Filt2Files), type="fastq")
Filt2 <- head(qaSummary[["readCounts"]])

dadaF1 <- readRDS(file.path(TnFpath1, "dadaF.rds"))
mergers1 <- readRDS(file.path(TnFpath1, "mergers.rds"))
seqtab1 <- readRDS(file.path(TnFpath1, "seqTab.rds"))
nochim1 <- readRDS(file.path(TnFpath1, "nochim.rds"))

sample.names <- sapply(strsplit(basename(files), "_"), `[`, 1)

print("Analysis alt 1")
getN <- function(x) sum(getUniques(x))
track <- cbind(NoAdap[, "read"], Trim[, "read"], Filt1[, "read"], Filt2[, "read"], sapply(dadaF1, getN), sapply(mergers1, getN), rowSums(seqtab1), rowSums(nochim1))
colnames(track) <- c("noAdap", "trimmed", "filtered", "matched", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
print(track)
