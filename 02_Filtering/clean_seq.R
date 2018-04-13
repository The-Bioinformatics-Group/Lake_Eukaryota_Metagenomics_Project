
## Removes all sequences in teh data containing any ambigous base (eg. and other letter than A, T, C, G)
## using the method clean
## Writes clean sequenses to files in new directory "cleaned"

#path <-  file.path("C:", "Users", "Karamech", "Documents", "Utbildning", "Kandidatarbete", "Data")
path <-  file.path("..", "00_Data") #Change to root of filtered Data
toClean <-  file.path(path, "trim_cutadapt") #Set to correct folder of files you want to clean
cleanPath <- file.path(path, "cleaned2") #path of new directory

library(ShortRead); packageVersion("ShortRead")

dir.create(cleanPath) #Creating directory for new files

##Finding all files and creating paths to new files with the same names
R1Files <- file.path(toClean, list.files(toClean, pattern="_R1_001.fastq.gz"))
R2Files <- file.path(toClean, list.files(toClean, pattern="_R2_001.fastq.gz"))

sample.names <- sapply(strsplit(basename(R1Files), "_"), `[`, 1)
R1CleanPath <- file.path(cleanPath, paste0(sample.names, "_R1_001.fastq.gz"))
R2CleanPath <- file.path(cleanPath, paste0(sample.names, "_R2_001.fastq.gz"))

##Looping over all R1 and R2 files respectively: read each file, clean and wite 
##to new fastq file 
for(i in 1:length(R1Files)){
	writeFastq(clean(readFastq(R1Files[i])), file=R1CleanPath[i], mode="w")
}
for(i in 1:length(R2Files)){
	writeFastq(clean(readFastq(R2Files[i])), file=R2CleanPath[i], mode="w")
}