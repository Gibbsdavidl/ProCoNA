\name{proconaPower}
\alias{proconaPower}
\title{proconaPower}
\description{this function returns the power used to build the procona network object}
\usage{proconaPower(p)}
\arguments{
  \item{p}{the procona obj}
}

\value{returns the power used to build the network}

\examples{
\dontrun{
data(ProCoNA_Data)
net1 <- buildProconaNetwork("peptide network", peptideData)
proconaPower(net1)
}
}

\author{David L Gibbs}





