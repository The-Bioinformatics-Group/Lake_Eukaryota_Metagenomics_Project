path <-  file.path("C:", "Users", "Karamech", "Documents", "Utbildning", "Kandidatarbete", "Data") #Change to root of all data files

library(dada2); packageVersion("dada2")

## Finding all data files in the root data directory
## and storing complete path in vector
dataFolders <- list.files(path) #vector of folders containing data
allFilesWithPath <- vector() #Initiating vector for all datafiles with complete paths
R1FilesWithPath <- vector()
R2FilesWithPath <- vector()  

findFilesInFolder <- function(folderpath){ #Finding data files in a folder/path (ening with .gz)
	gzfiles <- list.files(folderpath, pattern=".gz")
	return(gzfiles)
}

dataSets <- vector(length=0) 

for(i in dataFolders){ #Looping over folders to find datafiles
	filepath <- file.path(path,i)
	files <- findFilesInFolder(filepath)
	for(j in files){ #Looping over files and adding to vector
		dataSets <- c(dataSets, j)
		allFilesWithPath <- c(allFilesWithPath, file.path(path, i, j))
	}
}


## Sorting forward(R1) and reverse(R2) reads into two different vectors
R1pos <- grep("R1", allFilesWithPath)
R1FilesWithPath <- allFilesWithPath[R1pos]

R2pos <- grep("R2", allFilesWithPath)
R2FilesWithPath <- allFilesWithPath[R2pos]

# Can visualize quality profile of reads
# plotQualityProfile(R1FilesWithPath)
# plotQualityProfile(R2FilesWithPath)


##filtering and trimming all files and creating subdirectory with filtered files.
sample.names <- sapply(strsplit(basename(R1FilesWithPath), "_"), `[`, 1)

filt_path <- file.path(path, "filtered") # Path for placing filtered files
R1FilteredPath <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
R2FilteredPath <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(R1FilesWithPath, R1FilteredPath, R2FilesWithPath, R2FilteredPath, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE, otherwise to TRUE
head(out)