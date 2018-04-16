## Sample inference of alredy filtered and trimmed data
## Preparing by learning error rates and performing dereplication
## All sequences must not contain ambiguous bases and be longer than 5 bases.

#path <-  file.path("C:", "Users", "Karamech", "Documents", "Utbildning", "Kandidatarbete", "Data")
path <-  file.path("..", "00_Data") #Path to root of all data
TnFpath <- file.path(path, "clean_and_matched_2") #Path to directory of sequences that have been trimmed and filtered

library(dada2); packageVersion("dada2")

#finding files from path
R1files <- file.path(TnFpath, list.files(TnFpath, pattern="R1_001.fastq.gz"))
R2files <- file.path(TnFpath, list.files(TnFpath, pattern="R2_001.fastq.gz"))

##Learning error rates of the filtered files(Computationally intesive)
errF <- learnErrors(R1files, multithread=TRUE)
saveRDS(errF, file.path(TnFpath, "errF.rds"))
errR <- learnErrors(R2files, multithread=TRUE)
saveRDS(errR, file.path(TnFpath, "errR.rds"))

# Can plot error rates
# plotErrors(errF, nominalQ=TRUE)
# plotErrors(errR, nominalQ=TRUE)

##Dereplicate the filtered files (removing duplicates of sequences)
derepF <- derepFastq(R1files, verbose=TRUE) #Reading all files to memory? May have to modify befor running all data
saveRDS(derepF, file.path(TnFpath, "derepF.rds"))
derepR <- derepFastq(R2files, verbose=TRUE)
saveRDS(derepR, file.path(TnFpath, "derepR.rds"))

##Sample Inference
dadaF <- dada(derepF, err=errF, multithread=TRUE)
saveRDS(dadaF, file.path(TnFpath, "dadaF.rds"))
dadaR <- dada(derepR, err=errR, multithread=TRUE)
saveRDS(dadaR, file.path(TnFpath, "dadaR.rds"))

##Concatunate paired reads
#mergers <- mergePairs(dadaF, derepF, dadaR, derepR, justConcatenate=TRUE, verbose=TRUE) 
