## Prints the number of reads in each sample after each step in the data processing
## for data set 2/ trimming method 2
## Not fully implemented

path <-  file.path("..", "00_Data") #Path to root of all data
TnFpath2 <- file.path(path, "clean_and_matched_2") #Path to directory of sequences that have been trimmed and filtered

files <- list.files(TnFpath2, pattern="R1_001.fastq.gz")

library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")

#Raw

#Adapters removed
AdapFiles <- dir("../00_Data/AdaptersRemoved", "*R1_001.fastq.gz")
qaSummary <- qa(file.path(path, "AdaptersRemoved", AdapFiles), type="fastq")
NoAdap <- head(qaSummary[["readCounts"]])

#Trimmed 1
TrimFiles <- dir("../00_Data/trim_catadapt", "*R1_001.fastq.gz")
qaSummary <- qa(file.path(path, "trim_catadapt", TrimFiles), type="fastq")
Trim <- head(qaSummary[["readCounts"]])

#Filtered 1
FiltFiles <- dir("../00_Data/clean_and_matched_2", "*R1_001.fastq.gz")
qaSummary <- qa(file.path(path, "clean_and_matched_2", FiltFiles), type="fastq")
Filt <- head(qaSummary[["readCounts"]])

dadaF2 <- readRDS(file.path(TnFpath2, "dadaF.rds"))
mergers2 <- readRDS(file.path(TnFpath2, "mergers.rds"))
seqtab2 <- readRDS(file.path(TnFpath2, "seqTab.rds"))
nochim2 <- readRDS(file.path(TnFpath2, "nochim.rds"))

sample.names <- sapply(strsplit(basename(files), "_"), `[`, 1)

print("Analysis alt 2")
getN <- function(x) sum(getUniques(x))
track <- cbind(NoAdap[, "read"], Trim[, "read"], Filt[, "read"], sapply(dadaF2, getN), sapply(mergers2, getN), rowSums(seqtab2), rowSums(nochim2))
colnames(track) <- c("noAdap", "trimmed", "matched", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
print(track)