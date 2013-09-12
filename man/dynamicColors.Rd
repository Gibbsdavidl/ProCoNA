\name{dynamicColors}
\alias{dynamicColors}
\title{dynamicColors}
\description{Peptide module assignments before the merging step.}
\usage{dynamicColors(p)}
\arguments{
  \item{p}{The procona network object.}
}

\value{A numeric vector of module assignments, corresponding to peptides(proconaNet).}

\examples{
data(ProCoNA_Data)
net1 <- buildProconaNetwork("peptide network", peptideData, pow=12)
dynamicColors(net1)
}
\author{David L Gibbs}





