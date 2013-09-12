\name{networkType}
\alias{networkType}
\title{networkType}
\description{The networkType of the proconaNet, signed or unsigned.}
\usage{networkType(p)}
\arguments{
  \item{p}{the proconaNet obj}
}

\value{Returns the correlation type of the network (signed or unsigned).}


\examples{
\dontrun{
data(ProCoNA_Data)
net1 <- buildProconaNetwork("peptide network", peptideData)
networkType(net1)
}
}

\author{David L Gibbs}





