\name{compareNetworksWithFishersExactTest}
\alias{compareNetworksWithFishersExactTest}
\title{compareNetworksWithFishersExactTest}
\description{Fisher's exact test is used pairwise on modules to compare two networks. 
The arguments to Fisher's exact test are given below.

n == number of entities in the network

m == number of entities in intersection of two modules

d1 == number of entities in module A but not in module B

d2 == number of entities in module B but not in module A

2x2 matrix for the test is then:
m  d1
d2 n-d1-d2-m

}

\usage{compareNetworksWithFishersExactTest(peps1, peps2, colors1, colors2, 
    title = "", net1label = "", net2label = "")}

\arguments{
  \item{peps1}{Nodes in network 1, character vector}
  \item{peps2}{Nodes in network 2, character vector}
  \item{colors1}{modules for net 1}
  \item{colors2}{modules for net 2}
  \item{title}{Plot title}
  \item{net1label}{xlabel}
  \item{net2label}{ylabel}
}

\value{
	Returns fishers exact test -log pvalues and overlap matrix showing the number
	of shared members for each pair of modules.
}

\examples{
data(ProCoNA_Data)
net1 <- buildProconaNetwork("peptide network", peptideData, pow=12)
net2 <- buildProconaNetwork("peptide network", peptideData, pow=6)
compareNetworksWithFishersExactTest(peptides(net1), peptides(net2),
mergedColors(net1), mergedColors(net2), "network comparison", "net1", "net2")
}

\author{David L Gibbs}





