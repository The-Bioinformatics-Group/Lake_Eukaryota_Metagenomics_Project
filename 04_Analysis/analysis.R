## Creating a phyloseq object ps out of a dada2 sequence table and assigned taxonomy
## Extracts the sequence table and taxa table as RDS files from the specified directory

#path <-  file.path("..", "00_Data") 				#Path to root of all data
TnFpath <- file.path(path, "clean_and_matched_2") 	#Path to directory containing RDS files

path <- file.path("C:", "Users", "Karamech", "Documents", "Utbildning", "Kandidatarbete", "Data")
TnFpath <- file.path(path, "clean_and_matched_2")

##load packages
if (!require("vegan")) {
   install.packages("vegan", dependencies = TRUE)
   library(vegan)
   }
if (!require("ggplot2")) {
   install.packages("ggplot2", dependencies = TRUE)
   library(ggplot2)
   }
if (!require("ape")) {
   install.packages("ape", dependencies = TRUE)
   library(ape)
   }
if (!require("reshape2")) {
   install.packages("reshape2", dependencies = TRUE)
   library(reshape2)
   }
if (!require("phyloseq")) {
   source('http://bioconductor.org/biocLite.R')
	biocLite('phyloseq')
   library(phyloseq)
   }
if (!require("Hmisc")) {
   install.packages("Hmisc", dependencies = TRUE)
   library(Hmisc)
   }
if (!require("ggbiplot")) {
   library(devtools)
	install_github("ggbiplot", "vqv")
	library(ggbiplot)
   }

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
	CatchmentLakeRatio=metadata["CLratio"], SecchiDepth=metadata["Siktdjup"], Chlorophyll=metadata["Kfyll"], CaKMgNa=metadata["CaKMgNa"],
	pH=metadata["pH"], Conductivity=metadata["Kond_25"], TOC=metadata["TOC"], TotalP=metadata["TotP"],
	TotalN=metadata["TotN_TNb"], Temperature=metadata["Vattentemperatur"]))

#Remove outlier, very large catchment area Sample_ID=89
vardata1 <- vardata
vardata1["l100e89_R1_001.fastq.gz", 7:15] <- NA

#Normalize Land use types by catchment area
vardata1[, "Agriculturalkm"] <- vardata[,"Agriculturalkm"]/vardata[,"Catchment"]
vardata1[, "Wetlandkm"] <- vardata[,"Wetlandkm"]/vardata[,"Catchment"]
vardata1[, "Urbankm"] <- vardata[,"Urbankm"]/vardata[,"Catchment"]
vardata1[, "Grassbarekm"] <- vardata[,"Grassbarekm"]/vardata[,"Catchment"]
vardata1[, "Con_fores"] <- vardata[,"Con_fores"]/vardata[,"Catchment"]
vardata1[, "Dec_forest"] <- vardata[,"Dec_forest"]/vardata[,"Catchment"]
vardata1[, "Forest"] <- vardata[,"Forest"]/vardata[,"Catchment"]

#Standardize data
vardata2 <- decostand(vardata1, "stand", na.rm=TRUE)

#Calculation Shannon, chao1 and ACE and adding to data matrix

cnorm = round(norm[,colSums(norm)!=0])

eR5<-as.data.frame(t(estimateR(cnorm)))  #chao1 and ACE
Shan <- diversity(cnorm, index = "shannon", MARGIN = 1, base = exp(1))                # shannon
S <- specnumber(cnorm)
Shannon <- Shan/log(S)

a5<-cbind(Shannon, eR5$S.chao1, eR5$S.ACE)
#allData <- as.data.frame(merge(a5,vardata2, by.x = "row.names", by.y = "row.names"))
allData <- cbind(a5,vardata2)
 names(allData)[2] <- "Chao1"
 names(allData)[3] <- "ACE"

#Adding metadata to phyloseq
psS1 <- merge_phyloseq(psS, allData)

#cnorm = norm_short[,colSums(norm_short)!=0]

###########################Correlation matrix####################################
allData1 <- na.omit(allData)
test <- as.matrix(allData1)
test2 <- test[,2:ncol(test)]
DF = as.numeric(test2)
#norm_short <- cnorm[rownames(norm) %in% rownames(allData1),]
## allData1 <- decostand(allData, "stand", na.rm=TRUE) #Normalize after adding diversity?
res <- rcorr(allData1, type = "spearman")
# Extract the correlation coefficients and p-values with res$r & res$P
#Visualize p-values
symnum(res$P, abbr.colnames = TRUE)

#nona_res_P <- na.omit(res$P)
#sig_res_P <- nona_res_P[!(nona_res_P$value>0.05),]
#sig_res_r <- melted_res_r[rownames(melted_res_r) %in% rownames(sig_res_p),]
reorder_cormat <- function(cormat){
# Use correlation between variables as distance
dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
cormat <-cormat[hc$order, hc$order]
}
# Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
 }

P <- res$P
r <- res$r

P <- reorder_cormat(P)
P <- get_upper_tri(P)

r <- reorder_cormat(r)
r <- get_upper_tri(r)

#Plot correlation - heat map 
library(reshape2)
melted_res_P <- melt(P)
melted_res_r <- melt(r)

#Extract onlynon-NA, significant korrelations
melted_nona_res_P <- na.omit(melted_res_P)
melted_sig_res_p <- melted_nona_res_P[!(melted_nona_res_P$value>0.05),]

melted_sig_res_r <- melted_res_r[rownames(melted_res_r) %in% rownames(melted_sig_res_p),]

pdf("../00_Data/clean_and_matched_2/corrPlot3.pdf")
ggplot(data = melted_sig_res_r, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Spearman\nCorrelation") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(data = melted_sig_res_r, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#ggplot(data = sig_res_r, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

corrplot(res$r, p.mat = res$P, sig.level = .05, order="FPC")


########### RDA #########################
allData1 <- na.omit(allData)
norm_short <- norm[rownames(norm) %in% rownames(allData1),]
data.rda <- rda(norm_short ~ pH + TOC + TotP + TotN_TNb + Latitude + Agriculturalkm, data = allData1)

rdasc <- scores(data.rda, choices=c(1,2), display=c("sp","wa","cn", "bp"), scaling=2)
rdasites  <- data.frame(rdasc$sites[,1:2])      
rdabp  <- data.frame(rdasc$biplot[,1:2])

rda.plot <- ggplot(rdasites, aes(x=RDA1, y=RDA2, color=allData1$Shannon_Div)) + 
geom_point(size=4) +
geom_hline(yintercept=0, linetype="dotted") +
geom_vline(xintercept=0, linetype="dotted") +
coord_fixed()

rda.biplot.colour <- rda.plot + geom_segment(data=rdabp, aes(x=0, xend=RDA1, y=0, yend=RDA2), color="grey42", arrow=arrow(length=unit(0.01,"npc"))) + geom_text(data=rdabp, aes(x=RDA1,y=RDA2,label=rownames(rdabp),hjust=0.5*(1-sign(RDA1)),vjust=0.5*(1-sign(RDA2))), color="black", size=6) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(color="Latitude") + scale_colour_gradient(low="#56B1F7", high="#132B43") + scale_x_continuous(limits=c(-0.75,1.25))

#variation partitioning for chemistry vs latitude vs. catchment characteristics (only Agriculture left)
eukvar <- varpart(norm_short, ~ pH + TOC + TotP + TotN_TNb, ~Latitude, ~Agriculturalkm, data=allData1)

plot(eukvar)

#pdf("../00_Data/clean_and_matched_2/rdaPlot1.pdf")

############ Diversity plot
pdf("../00_Data/clean_and_matched_2/divPlot1.pdf")
plot(diversity$Latitude,diversity$Shannon_Div)
plot(diversity$Latitude,diversity$J5)
dev.off()

######################### NMDS ####################################
euk.nmds<-metaMDS(cnorm ,distance="bray", trymax=200)
pdf("../00_Data/clean_and_matched_2/NMDSplot1.pdf")
plot(euk.nmds, display = "sites")
dev.off()
n99<-scores(euk.nmds)
euk.nmds$stress
res1<-as.data.frame(merge(n99,allData1, by.x = "row.names", by.y = "row.names"))

pdf("../00_Data/clean_and_matched_2/AllNMDSplot2.pdf")
ggplot(res1, aes(NMDS2, NMDS1)) + geom_point(aes(size=2,colour=Shannon_Div)) + theme_bw(14) + xlab("Dimension 1") +ylab("Dimension 2")
dev.off()


######################################### PCA Analysis ##########################################
library(stats)
library(ggbiplot)
data <- allData1[, 1:23]
lakes <- allData1[, 1]
 
data.pca <- prcomp(data, center = TRUE, scale=TRUE) 

print(data.pca)

plot(data.pca, type = "l")

g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, groups = NULL, ellipse = TRUE, circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
---------------
theta <- seq(0,2*pi,length.out = 100)
circle <- data.frame(x = cos(theta), y = sin(theta))
p <- ggplot(circle,aes(x,y)) + geom_path()

loadings <- data.frame(data.pca$rotation, 
                       .names = row.names(data.pca$rotation))
p <- p + geom_text(data=loadings, 
              mapping=aes(x = PC1, y = PC2, label = .names, colour = .names)) +
  coord_fixed(ratio=1) +
  labs(x = "PC1", y = "PC2")


######################## Phyloseq tests ################################
#psS1.na = subset_taxa(psS1, is.na(Family))
#tax_table(psS1.na)[1:10, 1:6]
#psS1.notna = subset_taxa(psS1, !is.na(Family))

dino <- subset_taxa(psS, Phylum=="Dinoflagellata")
#plot_tree(dino, color="Order")

physeq = prune_taxa(taxa_names(psS)[1:200], psS)

top20 <- names(sort(taxa_sums(psS), decreasing=TRUE))[1:20]
psS.top20 <- transform_sample_counts(psS, function(OTU) OTU/sum(OTU))
psS.top20 <- prune_taxa(top20, psS.top20)

psS1.na = subset_taxa(psS.top20, is.na(Class))
#tax_table(psS1.na)[1:10, 1:6]
psS1.notna = subset_taxa(psS.top20, !is.na(Class))
plot_bar(psS1.na, psS1.notna, x="Phylum", color="Class")

#pdf("../00_Data/mothur_out_2/radTreePlot1.pdf")
#plot_tree(psS, nodelabf=nodeplotboot(60,60,3), ladderize="left") + coord_polar(theta="y")