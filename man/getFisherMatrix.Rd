\name{getFisherMatrix}
\alias{getFisherMatrix}
\title{getFisherMatrix}
\description{Fisher's exact test pairwise on modules.}
\usage{getFisherMatrix(net1, net2, colors1, colors2)}
\arguments{
  \item{net1}{Net1 list of entities in the network (nodes of network 1)}
  \item{net2}{Net2 list of entities in the network (nodes of network 2)}
  \item{colors1}{the module assignments for network 1}
  \item{colors2}{the module assignments for network 2}
}

\value{Returns the fisher test pvalues and count of overlapping peptides.}

\examples{
data(ProCoNA_Data)
net1 <- buildProconaNetwork("peptide network", peptideData, pow=12)
net2 <- buildProconaNetwork("peptide network", peptideData +
0.3*rnorm(length(peptideData)), pow=12)
getFisherMatrix(peptides(net1), peptides(net2), mergedColors(net1), mergedColors(net2))
}


\author{David L Gibbs}





