\name{mergedColors}
\alias{mergedColors}
\title{mergedColors}
\description{The vector of module assignments after the merging step.}
\usage{mergedColors(p)}
\arguments{
  \item{p}{The procona network object.}
}

\value{A numeric vector of module assignments corresponding to peptides(proconaNet).}

\examples{
data(ProCoNA_Data)
net1 <- buildProconaNetwork("peptide network", peptideData)
mergedColors(net1)
}

\author{David L Gibbs}





