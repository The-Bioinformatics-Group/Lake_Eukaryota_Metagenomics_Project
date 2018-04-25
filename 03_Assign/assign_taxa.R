## Assigning taxonomy to sequence table saved as rDS file in specified directory.
## Uses silva v132 training set, db located at root data path
## Saves taxa table and taxa.print table (without sequences in rownames) as RDS files
## in the same directory

path <-  file.path("..", "00_Data") #Path to root of all data
TnFpath <- file.path(path, "clean_and_matched_1") #Path to directory containing sequence table

##load package
library(dada2); packageVersion("dada2")

##Reading sequence table from file
nochim <- readRDS(file.path(TnFpath, "nochim.rds"))

##Assigning taxonomy
taxa <- assignTaxonomy(nochim, file.path(path, "silva_nr_v132_train_set.fa.gz"), multithread=TRUE)

##Make copy without sequences as rownames, print copy
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
print(taxa.print)

#Save both versions to RDS files
saveRDS(taxa, file.path(TnFpath, "taxa.rds"))
saveRDS(taxa.print, file.path(TnFpath, "taxa_print.rds"))
