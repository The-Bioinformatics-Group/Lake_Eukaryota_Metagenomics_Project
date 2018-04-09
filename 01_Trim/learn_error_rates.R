filtpath <-  file.path("..", "00_Data", "trimmed") #Change to root of filtered Data

library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")

R1Files <- file.path(filtpath, list.files(filtpath, pattern="_R1_"))
R2Files <- file.path(filtpath, list.files(filtpath, pattern="_R2_"))

sample.names <- sapply(strsplit(basename(R1Files), "_"), `[`, 1)
#filt_path <- file.path(path, "filtered") # Path for placing filtered files
R1CleanedPath <- file.path(path, "clean", paste0(sample.names, "_R1_clean.fastq.gz"))
R2CleanedPath <- file.path(path, "clean", paste0(sample.names, "_R2_clean.fastq.gz"))

#out <- filterAndTrim(R1Files, R1FilteredPath, R2Files, R2FilteredPath, truncLen=c(240,160),
 #             maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
  #            compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE, otherwise TRUE
#head(out)

for(i in 1:length(R1Files)){
	writeFastq(clean(readFastq(R1Files[i])), file=R1CleanedPath[i], mode="w")
}

for(i in 1:length(R2Files)){
	writeFastq(clean(readFastq(R2Files[i])), file=R2CleanedPath[i], mode="w")
}

cleanpath <- file.path(path, "clean")

R1Clean <- file.path(cleanpath, list.files(cleanpath, pattern="_R1_"))
R2Clean <- file.path(cleanpath, list.files(cleanpath, pattern="_R2_"))

##Learning error rates of the filtered files(Computationally intesive)
errF <- learnErrors(R1Clean, multithread=FALSE)
errR <- learnErrors(R2Clean, multithread=FALSE)

# Can plot error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)