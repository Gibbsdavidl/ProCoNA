\name{newProconaObj}
\alias{newProconaObj}
\title{newProconaObj}
\description{This function returns a peptide co-expression network object.}
\usage{
newProconaObj(networkName,
pepdat,
pow = NULL,
signed = F,
scaleFreeThreshold = 0.8,
deepSplit = 2,
minModuleSize = 30,
mergeThreshold = 0.1,
clusterType = "average",
performTOPermtest = TRUE,
toPermTestPermutes = 10000)
}
\arguments{
  \item{networkName}{Name of this network}
  \item{pepdat}{This variable is the data set with rows as samples and cols as peptides}
  \item{pow}{The scaling power, NULL if unknown}
  \item{signed}{Whether the sign is considered in constructing adjacency and TOM}
  \item{scaleFreeThreshold}{The threshold for fitting to scale-free topology.. will use closest power.}
  \item{deepSplit}{Course grain control of module size}
  \item{minModuleSize}{The minimum module size allowed}
  \item{mergeThreshold}{Below this threshold, modules are merged.}
  \item{clusterType}{Clustering option}
  \item{performTOPermtest}{Performs permutation testing on modules}
  \item{toPermTestPermutes}{Number of permutations to do.}
}

\value{returns the procona network object}

\author{David L Gibbs}





