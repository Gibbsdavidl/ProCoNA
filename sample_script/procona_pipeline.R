# Scripts by David L Gibbs
# gibbsd@ohsu.edu
# May 23, 2012

library(procona)
library(WGCNA)
library(reshape2)

options(stringsAsFactors=F)

# the mass tag db as a dataframe
masstagdb <- read.csv("MrOS_PeptideInfo_P724.csv")

#load phenotype
phenotype <- read.csv("mininet_pheno.csv")
names(phenotype)

#load peptides... and reshape them into a matrix...
# This is the point where samples should be removed for
# co-morbidities and etc..
peptideLongTable <- read.table("Mini_MS2_VP_Dave.txt",sep="\t",header=T)
peptideMatrix <- acast(data=peptideLongTable,
                       formula=Sample_ID~MassTagID,
                       value.var="Norm_UMCAbundance")


# subset the peptides to use only those with <= X% missing
pepDat <- log10(subsetPeptideData(peptideMatrix, percentageNAsAllowed=0.2))

# The order of the samples in the matrix should match the order
# of samples in the phenotype matrix.
print("MAKE SURE THE DATA ROWS MATCH THE PHENOTYPE MATRIX ROWS!!")

# to get a sense of what power would be best..
# 5 looks pretty good here.  Usually I take the first one over 0.8 or 0.85
pickSoftThreshold(pepDat, networkType="signed")

# make the procona object ... the network
# Try different settings of deepSplit (1-4), 2 is default, higher is more splitting 
miniNet <- newProconaObj("mini network", pepDat, signed="signed", pow=12, pearson=F)

# the significance of modules
miniNet@permtest

# the samples used.. in the order they came in
miniNet@samples[1:10]

# the peptides used...
miniNet@peptides[1:10]

# The module assignments to each peptide
miniNet@mergedColors[1:10]

# The module Summaries... eigenvectors for each module
miniNet@mergedMEs[1:5,1:5]

# Plot the dendrogram
plotDendroAndColors(miniNet@pepTree, miniNet@mergedColors, dendroLabels=F)

# correlation of module eigenvectors with phenotype
# returns data frame of correlations and pvalues
phenotypeCor <- correlationWithPhenotypesHeatMap(modules=1:5,   ## subset modules here
                                                 phenotypes=phenotype[,2:5], ## subset phentypes here
                                                 pdat=pepDat,
                                                 plotName=NULL, ## If you give a file name, plot is saves, not shown
                                                 pnet=miniNet,
                                                 textSize=0.7,
                                                 title="")
                 
# and getting specific module data out...

pepcor <- moduleMemberCorrelations(pnet=miniNet, pepdat=pepDat, phenotypes=phenotype)

#########################################################################
# quick function to write out the tables for specific modules.
moduleData <- function(pepnet, pepcors, module, pepinfo, fileprefix) {
  moduleX <- pepnet@peptides[which(pepnet@mergedColors==module)]
  moduleInfo <- pepinfo[which(pepinfo$Mass_Tag_ID %in% moduleX),]
  moduleCors <- pepcors[which(pepcors$Module==module),]
  corname <- paste(fileprefix, "_correlations.csv", sep="")
  write.table(moduleCors, file=corname, sep=",", row.names=F)
  infoname <- paste(fileprefix, "_peptide_info.csv", sep="")
  write.table(moduleInfo, file=infoname, sep=",", row.names=F)
}
########################################################################

moduleData(miniNet, phenotypeCor, 1, masstagdb, "Module_1")
