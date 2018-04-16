## Test for testing learnErrors
## and analysing the the error models through plotErrors
## All sequences must not contain ambiguous bases and be longer than 5 bases

#path <-  file.path("C:", "Users", "Karamech", "Documents", "Utbildning", "Kandidatarbete", "Data")
path <-  file.path("..", "00_Data") #Change to root of filtered Data
learnPath <- file.path(path, "cleaned2")

library(dada2); packageVersion("dada2")

R1learn <- file.path(learnPath, list.files(learnPath, pattern="R1_001.fastq.gz"))
R2learn <- file.path(learnPath, list.files(learnPath, pattern="R2_001.fastq.gz"))

##Learning error rates of the filtered files(Computationally intesive)
errF <- learnErrors(R1learn, multithread=TRUE)
errR <- learnErrors(R2learn, multithread=TRUE)

# Can plot error rates
pdf(file.path(learnPath, "errors.pdf"))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()
