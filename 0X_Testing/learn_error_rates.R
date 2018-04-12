#path <-  file.path("C:", "Users", "Karamech", "Documents", "Utbildning", "Kandidatarbete", "Data")
path <-  file.path("..", "00_Data") #Change to root of filtered Data
cleanpath <- file.path(path, "cleaned")

library(dada2); packageVersion("dada2")

R1Clean <- file.path(cleanpath, list.files(cleanpath, pattern="_R1_"))
R2Clean <- file.path(cleanpath, list.files(cleanpath, pattern="_R2_"))

##Learning error rates of the filtered files(Computationally intesive)
errF <- learnErrors(R1Clean, multithread=FALSE)
errR <- learnErrors(R2Clean, multithread=FALSE)

# Can plot error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)