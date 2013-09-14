\name{compareNetworksWithFishersExactTestProcona}
\alias{compareNetworksWithFishersExactTestProcona}
\alias{compareNetworksWithFishersExactTestProcona,proconaNet,proconaNet,character-method}
\title{compareNetworksWithFishersExactTestProcona}

\description{Convienence function for calling the compareNetworksWithFishersExactTest
using only two procona objects.}

\usage{compareNetworksWithFishersExactTestProcona(net1, net2, 
    title)}

\arguments{
  \item{net1}{procona object for network 1}
  \item{net2}{procona object for network 2}
  \item{title}{plot title}
}

\value{Returns a list of fisher -log pvalues, and overlaps between modules.}

\examples{
\dontrun{
data(ProCoNA_Data)
#net1 <- buildProconaNetwork("peptide network", peptideData)
#net2 <- buildProconaNetwork("peptide network", peptideData)
compareNetworksWithFishersExactTestProcona(net1, net2, "new comparison")
}
}

\author{David L Gibbs}





