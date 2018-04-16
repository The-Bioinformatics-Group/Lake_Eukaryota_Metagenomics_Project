path <-  file.path("C:", "Users", "Karamech", "Documents", "Utbildning", "Kandidatarbete", "Data")
path <-  file.path("..", "00_Data") #Change to root of filtered Data
TrimPath <- file.path(path, "trimmedTailsw")

library(dada2); packageVersion("dada2")

R1files <- file.path(TrimPath, list.files(TrimPath, pattern="_R1_001.fastq.gz"))
R2files <- file.path(TrimPath, list.files(TrimPath, pattern="_R2_001.fastq.gz"))
sample.names <- sapply(strsplit(basename(R1files), "_"), `[`, 1)

dir.create(file.path(path, "clean_and_matched"))

R1newPath <- file.path(path, "clean_and_matched", paste0(sample.names, "_R1_001.fastq.gz"))
R2newPath <- file.path(path, "clean_and_matched", paste0(sample.names, "_R2_001.fastq.gz"))

for(i in 1:length(R1files)){
	pair <- c(R1files[i], R2files[i])
	out <- c(R1newPath[i], R2newPath[i])
	fastqPairedFilter(pair, out, matchIDs=TRUE)
	print(out)
}

