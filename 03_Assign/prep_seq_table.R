## Makes a sequnce table from dada2 mergers object and removes chimeras
## Reads merger from RDS file located in speciefied directory
## Sequence table and Sequence table with chimeras removed are saved as
## separate RDS files to the same directory

path <-  file.path("..", "00_Data") #Path to root of all data
TnFpath <- file.path(path, "clean_and_matched_2") #Path to directory containing RDS file with mergers object

##Load package
library(dada2); packageVersion("dada2")

##Read mergers object from file
mergers <- readRDS(file.path(TnFpath, "mergers.rds"))

##Make sequence table
seqtab <- makeSequenceTable(mergers)

##Remove bimeras from sequence table and save as new table
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#dim(seqtab.nochim)

##Save both sequence tables to file
saveRDS(seqtab, file.path(TnFpath, "seqTab.rds"))
saveRDS(seqtab.nochim, file.path(TnFpath, "nochim.rds"))
