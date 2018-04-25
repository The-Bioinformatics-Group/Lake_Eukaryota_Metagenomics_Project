## Trimming poor quality data using trimTails. 
## Uses sliding window 2*halfwith+1, if a window contains k bases with phred 
## score <= a the sequence is cut at the start of the window
## Uses the dataset where primers have already been removed, creates new directory for trimmed files
## One of two method used for trimming the data. Produces sequences of variable length

#path <-  file.path("C:", "Users", "Karamech", "Documents", "Utbildning", "Kandidatarbete", "Data")
path <-  file.path("..", "00_Data") #Change to root of filtered Data
NoAdapPath <-  file.path(path, "AdaptersRemoved") #Set to directory whith files to trim
TrimPath <- file.path(path, "trimmedTailsw") #new directory with trimmed sequences, created in root for data

library(ShortRead); packageVersion("ShortRead")

dir.create(TrimPath) #Creating directory for new files

##Finding all files and creating paths to new files with the same names
R1Files <- file.path(NoAdapPath, list.files(NoAdapPath, pattern="_R1_001.fastq.gz"))
R2Files <- file.path(NoAdapPath, list.files(NoAdapPath, pattern="_R2_001.fastq.gz"))

sample.names <- sapply(strsplit(basename(R1Files), "_"), `[`, 1)
R1TrimPath <- file.path(TrimPath, paste0(sample.names, "_R1_001.fastq.gz"))
R2TrimPath <- file.path(TrimPath, paste0(sample.names, "_R2_001.fastq.gz"))

## Filtering files using trimTailsw.
## Same quality check for R1 and R2
trimWhen <- character(length=1)
trimWhen[1] <- "+" # ASCII character, trimTailsw trims at first instance from the left of base with phred score <= trimWhen"
nbrFail <- integer(1)
nbrFail[1] <- 5 # Number of bases with phred score <= to trimWhen nessecary for triggering trimming
trimTailw(R1Files, k=nbrFail, a=trimWhen, halfwidth=5, destinations=R1TrimPath, ranges=FALSE, right=TRUE)
trimTailw(R2Files, k=nbrFail, a=trimWhen, halfwidth=5, destinations=R2TrimPath, ranges=FALSE, right=TRUE)
