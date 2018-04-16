## Sample inference of alredy filtered and trimmed data
## Preparing by learning error rates and performing dereplication
## All sequences must not contain ambiguous bases and be longer than 5 bases.

#path <-  file.path("C:", "Users", "Karamech", "Documents", "Utbildning", "Kandidatarbete", "Data")
path <-  file.path("..", "00_Data") #Path to root of all data
TnFpath <- file.path(path, "cleaned2") #Path to directory of sequences that have been trimmed and filtered

library(dada2); packageVersion("dada2")

#finding files from path
R1files <- file.path(TnFpath, list.files(TnFpath, pattern="R1_001.fastq.gz"))
R2files <- file.path(TnFpath, list.files(TnFpath, pattern="R2_001.fastq.gz"))

##Learning error rates of the filtered files(Computationally intesive)
errF <- learnErrors(R1files, multithread=FALSE)
errR <- learnErrors(R2files, multithread=FALSE)

# Can plot error rates
# plotErrors(errF, nominalQ=TRUE)
# plotErrors(errR, nominalQ=TRUE)

##Dereplicate the filtered files (removing duplicates of sequences)
derepF <- derepFastq(R1files, verbose=TRUE) #Reading all files to memory? May have to modify befor running all data
derepR <- derepFastq(R2files, verbose=TRUE)

names(derepF) <- sample.names
names(derepR) <- sample.names

##Sample Inference
dadaF <- dada(derepF, err=errF, multithread=FALSE)
dadaR <- dada(derepR, err=errR, multithread=FALSE)

##Concatunate paired reads

mergers <- mergePairs(dadaF, derepF, dadaR, derepR, justConcatenate=TRUE, verbose=TRUE) 
#Should merge F and R reads but does not, need to check if read overlap default for function is min 20 nts