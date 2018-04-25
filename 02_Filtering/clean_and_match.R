## Removes all sequences containing any ambiguos base (N).
## Tries to match up all to pairs of forward and reverse reads by ID,  
## removing all sequences that have been orphaned by filtering.
## Creates new directory for clean and matched files 

path <-  file.path("..", "00_Data") #Change to root of filtered Data
TrimPath <- file.path(path, "longSeq") #Directory of trimmed (and filtered) files

##Load library
library(dada2); packageVersion("dada2")

##Finding all files in directory
R1files <- file.path(TrimPath, list.files(TrimPath, pattern="_R1_001.fastq.gz"))
R2files <- file.path(TrimPath, list.files(TrimPath, pattern="_R2_001.fastq.gz"))
sample.names <- sapply(strsplit(basename(R1files), "_"), `[`, 1)

##Create new directory
dir.create(file.path(path, "clean_and_matched_1"))

##Path and names of new fles
R1newPath <- file.path(path, "clean_and_matched_1", paste0(sample.names, "_R1_001.fastq.gz"))
R2newPath <- file.path(path, "clean_and_matched_1", paste0(sample.names, "_R2_001.fastq.gz"))

##For every pair of files, do paired filtering and output to new pair of files
for(i in 1:length(R1files)){
	pair <- c(R1files[i], R2files[i])
	out <- c(R1newPath[i], R2newPath[i])
	fastqPairedFilter(pair, out, matchIDs=TRUE) #matchIDs=TRUE removes seq with no matching ID in opposite direction,
	print(out)									#max allowed N defaults to 0
}
