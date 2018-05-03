##

path <-  file.path("..", "00_Data") 				#Path to root of all data
TnFpath <- file.path(path, "clean_and_matched_1") 	#Path to directory containing RDS files

library(ShortRead); packageVersion("ShortRead")
library(dada2); packageVersion("dada2")

nochim <- readRDS(file.path(TnFpath, "nochim.rds"))
seq <- colnames(nochim)

prefix <- "seq"
suffix <- 1:length(seq)
names <- paste0(prefix, suffix)
names(seq) <- names

writeFasta(seq, file.path(TnFpath, "seq.fasta"))