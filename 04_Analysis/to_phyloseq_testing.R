## Creating a phyloseq object ps out of a dada2 sequence table and assigned taxonomy
## Extracts the sequence table and taxa table as RDS files from the specified directory

path <-  file.path("..", "00_Data") 				#Path to root of all data
TnFpath <- file.path(path, "clean_and_matched_2_old") 	#Path to directory containing RDS files

##load packages
library(phyloseq); packageVersion("phyloseq")
#library(ggplot2); packageVersion("ggplot2")

##Reading sequence table and taxa table
nochim <- readRDS(file.path(TnFpath, "nochim.rds"))
taxa <- readRDS(file.path(TnFpath, "taxa.rds"))

##Formating sample data. Metadata can/should be included here.
samples.out <- rownames(nochim)
subject <- sapply(strsplit(samples.out, "_"), `[`, 1)
samdf <- data.frame(Subject=subject)
rownames(samdf) <- samples.out

##Creating and saving phyloseq object
ps <- phyloseq(otu_table(nochim, taxa_are_rows=FALSE), sample_data(samdf), tax_table(taxa))
#saveRDS(ps, file.path(TnFpath, "ps.rds"))

metafile <- read.csv(file.path(path, "100lakes_metadata_ed_short.csv"), sep=";")
metadata <- array(metafile, dimnames=list(sample_names(ps), colnames(metafile)))

vardata = sample_data(data.frame(Sample=metadata["sampleID"], Latitude=metadata["Latitude"], pH=metadata["pH"]))
ps1 <- merge_phyloseq(ps, vardata)

psS1.na = subset_taxa(ps1, is.na(Family))
#tax_table(psS1.na)[1:10, 1:6]
psS1.notna = subset_taxa(ps1, !is.na(Family))