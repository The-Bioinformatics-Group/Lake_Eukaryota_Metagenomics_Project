path <-  file.path("..", "00_Data") #Change to root of all data files

library(ShortRead); packageVersion("ShortRead")

## Finding all data files in the root data directory
## and storing complete path in vector
dataFolders <- list.files(path, pattern="Sample") #vector of folders containing data
allFilesWithPath <- vector() #Initiating vector for all datafiles with complete paths
R1FilesWithPath <- vector()
R2FilesWithPath <- vector()  

findFilesInFolder <- function(folderpath){ #Finding data files in a folder/path (ening with .gz)
	gzfiles <- list.files(folderpath, pattern=".gz")
	return(gzfiles)
}

dataSets <- vector(length=0) 

stop <- 0
for(i in dataFolders){ #Looping over folders to find datafiles
	stop <- stop+1
	if(stop > 5){ #Change according to number of samples you want to run
		break
	}
	filepath <- file.path(path,i)
	files <- findFilesInFolder(filepath)
	for(j in files){ #Looping over files and adding to vector
		dataSets <- c(dataSets, j)
		allFilesWithPath <- c(allFilesWithPath, file.path(path, i, j))
	} 
}

print(allFilesWithPath)


## Sorting forward(R1) and reverse(R2) reads into two different vectors
R1pos <- grep("R1", allFilesWithPath)
R1FilesWithPath <- allFilesWithPath[R1pos]

R2pos <- grep("R2", allFilesWithPath)
R2FilesWithPath <- allFilesWithPath[R2pos]

## Create new directory and names for filtered files
sample.names <- sapply(strsplit(basename(R1FilesWithPath), "_"), `[`, 1)
trim_path <- file.path(path, "trimmed") #Path for placing copies to be trimmed files
dir.create(trim_path)
R1TrimPath <- file.path(trim_path, paste0(sample.names, "_R1_trim.fastq.gz"))
R2TrimPath <- file.path(trim_path, paste0(sample.names, "_R2_trim.fastq.gz"))

## Filtering files using trimTailsw
trimWhen <- character(length=1)
trimWhen[1] <- "0" # ASCII character, trimTailsw trims at first instance from the left of base with phred score <= trimWhen"
nbrFail <- integer(1)
nbrFail[1] <- 3 # Number of bases with phred score <= to trimWhen nessecary for triggering trimming
trimTailw(R1FilesWithPath, k=nbrFail, a=trimWhen, halfwidth=5, destinations=R1TrimPath, ranges=FALSE, right=TRUE)
trimTailw(R2FilesWithPath, k=nbrFail, a=trimWhen, halfwidth=5, destinations=R2TrimPath, ranges=FALSE, right=TRUE)