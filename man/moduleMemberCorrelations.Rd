\name{moduleMemberCorrelations}
\alias{moduleMemberCorrelations}
\title{moduleMemberCorrelations}

\description{Computes the relation between peptides and eigenvector summaries
and also peptides and phenotypes.}

\usage{moduleMemberCorrelations(pnet, pepdat, phenotypes)}

\arguments{
  \item{pnet}{The peptide net object}
  \item{pepdat}{The peptide data matrix}
  \item{phenotypes}{The matrix of traits}
}

\value{Matrix of Pearson correlations with peptides in rows.}

\examples{
data(ProCoNA_Data)
net1 <- buildProconaNetwork("peptide network", peptideData)
n <- length(samples(net1))
phenotypes <- matrix(rnorm(10*n), nrow=60)
pepcor <- moduleMemberCorrelations(net1, peptideMatrix, phenotypes)

# To plot the heatmap: 
# moduleCors <- correlationWithPhenotypesHeatMap(net1, phenotypes, modules = 1:5, 
#    plot = NULL, title = "Module-trait relationships", textSize = 0.5)

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

# WRITE OUT A TABLE WITH THE BELOW FUNCTION CALL  :)#
# moduleData(peptideNetwork, pepcor, 1, masstagdb, "Module_1")
}
\author{David L Gibbs}





