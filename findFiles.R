
path <-  file.path("C:", "Users", "Karamech", "Documents", "Utbildning", "Kandidatarbete", "Data") 

library(dada2); packageVersion("dada2")

dataFolders <- list.files(path) 
allFilesWithPath <- vector() 

findFilesInFolder <- function(folderpath){ 
	gzfiles <- list.files(folderpath, pattern=".gz")
	return(gzfiles)
}

dataSets <- vector(length=0) 

for(i in dataFolders){ 
	filepath <- file.path(path,i)
	files <- findFilesInFolder(filepath)
	for(j in files){
		dataSets <- c(dataSets, j)
		allFilesWithPath <- c(allFilesWithPath, file.path(path, i, j))
	}

}

R1FilesWithPath <- vector()
R2FilesWithPath <- vector()

for (i in allFilesWithPath) {
	R1FilesWithPath <- 
}
