## Concatenates pairs of reverse and forward read at the 3' end
## inserting 10xN at the concatenation site
## Uses dada2 objects saved as RDS objects in specified directory
## Saves the merging to RDS file in the same directory

path <-  file.path("..", "00_Data") #Path to root of all data
TnFpath <- file.path(path, "clean_and_matched_2") #Path to directory containing dada2 objects as RDS files

##Load package
library(dada2); packageVersion("dada2")

##Read dada2 objects from file
dadaF <- readRDS(file.path(TnFpath, "dadaF.rds"))
dadaR <- readRDS(file.path(TnFpath, "dadaR.rds"))
derepF <- readRDS(file.path(TnFpath, "derepF.rds"))
derepR <- readRDS(file.path(TnFpath, "derepR.rds"))

##Concatenate
mergers <- mergePairs(dadaF, derepF, dadaR, derepR, justConcatenate=TRUE, verbose=TRUE)

##Save mergers to RDS file
saveRDS(mergers, file.path(TnFpath, "mergers.rds"))

