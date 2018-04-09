filtpath <-  file.path(file.path("..", "00_Data", "trimmed") #Change to root of filtered Data

library(dada2); packageVersion("dada2")

R1Files <- file.path(filtpath, list.files(filtpath, pattern="_R1_"))
R2Files <- file.path(filtpath, list.files(filtpath, pattern="_R2_"))

##Learning error rates of the filtered files(Computationally intesive)
errF <- learnErrors(R1Files, multithread=TRUE)
errR <- learnErrors(R2Files, multithread=TRUE)

# Can plot error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)