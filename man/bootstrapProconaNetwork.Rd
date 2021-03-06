\name{bootstrapProconaNetwork}
\alias{bootstrapProconaNetwork}
\title{bootstrapProconaNetwork}

\description{This function returns a peptide co-expression network object based on a bootstrapped correlation matrix.}

\usage{bootstrapProconaNetwork(networkName = "bootstrap procona", pepdat = NULL, 
    pow = NULL, powMax = 20, networkType = "signed", scaleFreeThreshold = 0.8, 
    deepSplit = 2, minModuleSize = 30, mergeThreshold = 0.1, 
    clusterType = "average", pamRespectsDendro = T, performTOPermtest = TRUE, 
    toPermTestPermutes = 100, bootstrapThreshold = 1e-04)}

\arguments{
  \item{networkName}{Name of this network}
  \item{pepdat}{This variable is the data set with rows as samples and cols as peptides}
  \item{pow}{The scaling power, NULL if unknown}
  \item{powMax}{The maximum power to be searched.}
  \item{networkType}{Whether the sign is considered in constructing adjacency and TOM}
  \item{scaleFreeThreshold}{The threshold for fitting to scale-free topology.. will use closest power.}
  \item{deepSplit}{Course grain control of module size}
  \item{minModuleSize}{The minimum module size allowed}
  \item{mergeThreshold}{Below this threshold, modules are merged.}
  \item{clusterType}{Clustering option}
  \item{pamRespectsDendro}{When cutting the dendrogram, pay attention to branch membership.}
  \item{performTOPermtest}{Performs permutation testing on modules}
  \item{toPermTestPermutes}{Number of permutations to do.}
  \item{bootstrapThreshold}{When to stop resampling...}
}

\value{returns the procona network object}

\examples{

data(ProCoNA_Data)
net <- bootstrapProconaNetwork("peptide network", peptideData,
performTOPermtest=FALSE, bootstrapThreshold=0.1)
}
\author{David L Gibbs}





