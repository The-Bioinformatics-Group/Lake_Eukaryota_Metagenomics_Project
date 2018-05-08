## Creating a phyloseq object ps out of a dada2 sequence table and assigned taxonomy
## Extracts the sequence table and taxa table as RDS files from the specified directory

path <-  file.path("..", "00_Data") 				#Path to root of all data
TnFpath <- file.path(path, "clean_and_matched_2") 	#Path to directory containing RDS files

##load packages
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(ape); packageVersion("ape")
library(ggplot2); packageVersion("ggplot2")
library("Hmisc")

##Reading sequence table and taxa table
nochim <- readRDS(file.path(TnFpath, "nochim.rds"))
taxa <- readRDS(file.path(TnFpath, "taxa.rds"))

##Normalize data
prob <- drarefy(nochim, sample=5666)
norm <- nochim*prob

##Formating sample data. Metadata can/should be included here.
samples.out <- rownames(norm)
subject <- sapply(strsplit(samples.out, "_"), `[`, 1)
samdf <- data.frame(Subject=subject)
rownames(samdf) <- samples.out

#Reading tree generated with fastTree
phylo <- read.tree(file.path(path, "mothur_out_2", "phylo.tree"))
#Fixing tip labels to match otu names
phylo$tip.label <- sapply(strsplit(phylo$tip.label, "q"), `[`, 2)
phylo$tip.label <- as.numeric(phylo$tip.label)
newlabels <- c()
for (i in phylo$tip.label){
	newlabels[i] <- rownames(taxa)[phylo$tip.label[i]]
}
phylo$tip.label <- newlabels

##Creating and saving phyloseq object
psS <- phyloseq(otu_table(norm, taxa_are_rows=FALSE), sample_data(samdf), tax_table(taxa), phy_tree(phylo))
#saveRDS(ps, file.path(TnFpath, "ps.rds"))


######################### Metadata and diversity ########################################################
metafile <- read.csv(file.path(path, "100lakes_metadata_ed_1.csv"), sep=";")
metadata <- array(metafile, dimnames=list(sample_names(psS), colnames(metafile)))

vardata = sample_data(data.frame(Latitude=metadata["Latitude"], Longitude=metadata["Longitude"], Aluminum=metadata["AlICPAES"], 
	AlkAcid=metadata["AlkAcid"], Turbidity=metadata["Turb_FNU"], Area=metadata["Area"], CatchmentArea=metadata["Catchment"], 
	Agriculture=metadata["Agriculturalkm"], Wetland=metadata["Wetlandkm"], Urban=metadata["Urbankm"], Grassland=metadata["Grassbarekm"], 
	ConiferousForest=metadata["Con_fores"], DeciduousForest=metadata["Dec_forest"], MixedForest=metadata["Forest"],
	CatchmentLakeRatio=metadata["CLratio"], SecchiDepth=metadata["Siktdjup"], Oxygen=metadata["Syrgashalt"], 
	pH=metadata["pH"], Conductivity=metadata["Kond_25"], TOC=metadata["TOC"], TotalP=metadata["TotP"],
	TotalN=metadata["TotN_TNb"], Temperature=metadata["Vattentemperatur"]))

vardata1 <- decostand(vardata, "stand", na.rm=TRUE)

#Shannon diversity
ShanDiv <- diversity(norm, index = "shannon", MARGIN = 1, base = exp(1))
MShanDiv <- matrix(ShanDiv)
rownames(MShanDiv) <- names(ShanDiv)
colnames(MShanDiv) <- "Shannon_Div"
allData = cbind(vardata1, MShanDiv)

psS1 <- merge_phyloseq(psS, allData)

######## Phyloseq tests ########################
#psS1.na = subset_taxa(psS1, is.na(Family))
#tax_table(psS1.na)[1:10, 1:6]
#psS1.notna = subset_taxa(psS1, !is.na(Family))

dino <- subset_taxa(psS, Phylum=="Dinoflagellata")
plot_tree(dino, color="Order")

physeq = prune_taxa(taxa_names(psS)[1:200], psS)

top20 <- names(sort(taxa_sums(psS), decreasing=TRUE))[1:50]
psS.top20 <- transform_sample_counts(psS, function(OTU) OTU/sum(OTU))
psS.top20 <- prune_taxa(top20, psS.top20)


##########################Correlation matrix####################################
## allData1 <- decostand(allData, "stand", na.rm=TRUE) #Normalize after adding diversity?
res <- rcorr(as.matrix(allData1), type = "spearman")
# Extract the correlation coefficients and p-values with res$r & res$P
#Visualize p-values
symnum(res$P, abbr.colnames = TRUE)

####### Plot correlation - heat map #######
library(reshape2)
melted_res <- melt(res$P)
ggplot(data = melted_res, aes(x=Var1, y=Var2, fill=value)) + geom_tile()