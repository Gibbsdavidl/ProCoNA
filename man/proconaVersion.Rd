\name{proconaVersion}
\alias{proconaVersion}
\title{proconaVersion}
\description{Keep track of the version of procona used to build networks.}
\usage{proconaVersion(p)}
\arguments{
  \item{p}{The procona network object.}
}
\examples{
\dontrun{
data(ProCoNA_Data)
net1 <- buildProconaNetwork("peptide network", peptideData)
proconaVersion(net1)
}
}

\author{David L Gibbs}





