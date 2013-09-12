\name{utri}
\alias{utri}
\title{utri}
\description{The upper triangle of a matrix}
\usage{utri(mat)}
\arguments{
  \item{mat}{A matrix}
}

\value{Returns a vector}
\examples{
data(ProCoNA_Data)
net1 <- buildProconaNetwork("peptide network", peptideData)
utri(TOM(net1))
}

\author{David L Gibbs}





