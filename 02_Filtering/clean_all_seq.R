
## Removes all sequences in teh data containing any ambigous base (eg. and other letter than A, T, C, G)
## using the method clean
## Writes clean sequenses to files in new directory "cleaned"

#path <-  file.path("C:", "Users", "Karamech", "Documents", "Utbildning", "Kandidatarbete", "Data")
path <-  file.path("..", "00_Data") #Change to root of filtered Data
toClean1 <-  file.path(path, "trimmedTailsw") #Set to correct folder of files you want to clean
toClean2 <-  file.path(path, "trim_cutadapt") #Set to correct folder of files you want to clean
cleanPath1 <- file.path(path, "cleaned1") #path of new directory
cleanPath2 <- file.path(path, "cleaned2") #path of new directory

library(ShortRead); packageVersion("ShortRead")

dir.create(cleanPath) #Creating directory for new files

##Finding all files and creating paths to new files with the same names
R1Files1 <- file.path(toClean1, list.files(toClean1, pattern="_R1_001.fastq.gz"))
R2Files1 <- file.path(toClean1, list.files(toClean1, pattern="_R2_001.fastq.gz"))
R1Files2 <- file.path(toClean2, list.files(toClean2, pattern="_R1_001.fastq.gz"))
R2Files2 <- file.path(toClean2, list.files(toClean2, pattern="_R2_001.fastq.gz"))

sample.names <- sapply(strsplit(basename(R1Files), "_"), `[`, 1)
R1CleanPath1 <- file.path(cleanPath1, paste0(sample.names, "_R1_001.fastq.gz"))
R2CleanPath1 <- file.path(cleanPath1, paste0(sample.names, "_R2_001.fastq.gz"))
R1CleanPath2 <- file.path(cleanPath2, paste0(sample.names, "_R1_001.fastq.gz"))
R2CleanPath2 <- file.path(cleanPath2, paste0(sample.names, "_R2_001.fastq.gz"))

##Looping over all R1 and R2 files respectively: read each file, clean and wite 
##to new fastq file 
for(i in 1:length(R1Files1)){
	writeFastq(clean(readFastq(R1Files1[i])), file=R1CleanPath1[i], mode="w")
	writeFastq(clean(readFastq(R1Files2[i])), file=R1CleanPath2[i], mode="w")
	writeFastq(clean(readFastq(R2Files1[i])), file=R2CleanPath1[i], mode="w")
	writeFastq(clean(readFastq(R2Files2[i])), file=R2CleanPath2[i], mode="w")
}