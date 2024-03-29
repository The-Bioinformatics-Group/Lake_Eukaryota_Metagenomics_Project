
path <-  file.path("..", "00_Data") #Change to root of all data files

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