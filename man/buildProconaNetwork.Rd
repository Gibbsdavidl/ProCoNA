\name{buildProconaNetwork}
\alias{buildProconaNetwork}
\title{buildProconaNetwork}
\description{This function returns a peptide co-expression network object.}
\usage{buildProconaNetwork(networkName = "ProCoNA", pepdat = NULL, pow = NULL, 
    powMax = 20, networkType = "signed", pearson = F, scaleFreeThreshold = 0.8, 
    deepSplit = 2, minModuleSize = 30, mergeThreshold = 0.1, 
    clusterType = "average", pamRespectsDendro = T, performTOPermtest = TRUE, 
    toPermTestPermutes = 100)}
\arguments{
  \item{networkName}{Name of this network}
  \item{pepdat}{This variable is the data set with rows as samples and cols as peptides}
  \item{pow}{The scaling power, NULL if unknown}
  \item{powMax}{The maximum power to be searched.}
  \item{networkType}{Should the sign be considered in constructing adjacency and TOM ("signed" or "unsigned")}
  \item{pearson}{use Pearson's cor or the robust bi-weight correlation}
  \item{scaleFreeThreshold}{The threshold for fitting to scale-free topology.. will use closest power.}
  \item{deepSplit}{Course grain control of module size}
  \item{minModuleSize}{The minimum module size allowed}
  \item{mergeThreshold}{Below this threshold, modules are merged.}
  \item{clusterType}{Clustering option}
  \item{pamRespectsDendro}{When cutting the dendrogram, pay attention to branch membership.}
  \item{performTOPermtest}{Performs permutation testing on modules}
  \item{toPermTestPermutes}{Number of permutations to do.}
}

\value{returns the procona network object}

\details{
The procona network object contains a number of slots 
which store information relevant to the construction of 
the network. Accessor functions provide direct access to 
the slots. See getSlots("proconaNet") for a complete list.
}

\examples{
data(ProCoNA_Data)
net <- buildProconaNetwork("peptide network", peptideData)
}
\author{David L Gibbs}





