\name{samples}
\alias{samples}
\title{samples}
\description{The sample names used to build the network. Useful for checking phenotype ordering.}
\usage{samples(p)}
\arguments{
  \item{p}{the procona obj}
}

\value{Returns the character vector of sample IDs.}
\examples{
\dontrun{
data(ProCoNA_Data)
net1 <- buildProconaNetwork("peptide network", peptideData)
samples(net1)
}
}
\author{David L Gibbs}





