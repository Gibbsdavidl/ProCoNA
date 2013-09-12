\name{TOM}
\alias{TOM}
\title{TOM}
\description{The topological overlap matrix.}
\usage{TOM(p)}
\arguments{
  \item{p}{the procona obj}
}

\value{Returns the topological overlap matrix.}
\examples{
data(ProCoNA_Data)
net1 <- buildProconaNetwork("peptide network", peptideData)
tomMat <- TOM(net1)
}

\author{David L Gibbs}





