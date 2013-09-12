\name{peptides}
\alias{peptides}
\title{peptides}
\description{The peptides used to build the peptideNet object.}
\usage{peptides(p)}
\arguments{
  \item{p}{The proconaNet object}
}

\value{Returns a vector of peptide names or IDs.}

\examples{
\dontrun{
data(ProCoNA_Data)
net1 <- buildProconaNetwork("peptide network", peptideData)
peptides(net1)
}
}

\author{David L Gibbs}





