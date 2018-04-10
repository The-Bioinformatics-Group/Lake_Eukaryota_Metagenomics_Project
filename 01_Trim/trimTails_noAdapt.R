#path <-  file.path("C:", "Users", "Karamech", "Documents", "Utbildning", "Kandidatarbete", "Data")
path <-  file.path("..", "00_Data") #Change to root of filtered Data
NoAdapPath <-  file.path(path, "AdaptersRemoved")
TrimPath <- file.path(path, "trimmedTailsw")
dir.create(TrimPath)

library(ShortRead); packageVersion("ShortRead")

R1Files <- file.path(NoAdapPath, list.files(NoAdapPath, pattern="_R1_"))
R2Files <- file.path(NoAdapPath, list.files(NoAdapPath, pattern="_R2_"))

sample.names <- sapply(strsplit(basename(R1Files), "_"), `[`, 1)
#filt_path <- file.path(path, "filtered") # Path for placing filtered files
R1TrimPath <- file.path(TrimPath, paste0(sample.names, "_R1_trimTails.fastq.gz"))
R2TrimPath <- file.path(TrimPath, paste0(sample.names, "_R2_trimTails.fastq.gz"))

## Filtering files using trimTailsw
trimWhen <- character(length=1)
trimWhen[1] <- "+" # ASCII character, trimTailsw trims at first instance from the left of base with phred score <= trimWhen"
nbrFail <- integer(1)
nbrFail[1] <- 5 # Number of bases with phred score <= to trimWhen nessecary for triggering trimming
trimTailw(R1Files, k=nbrFail, a=trimWhen, halfwidth=5, destinations=R1TrimPath, ranges=FALSE, right=TRUE)
trimTailw(R2Files, k=nbrFail, a=trimWhen, halfwidth=5, destinations=R2TrimPath, ranges=FALSE, right=TRUE)