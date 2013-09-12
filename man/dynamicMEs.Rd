\name{dynamicMEs}
\alias{dynamicMEs}
\title{dynamicMEs}
\description{The module eigenvectors before the module merging step.}
\usage{dynamicMEs(p)}
\arguments{
  \item{p}{The procona network object.}
}

\value{Returns the matrix of module eigenvectors, before module merging, in columns.}

\examples{
data(ProCoNA_Data)
net1 <- buildProconaNetwork("peptide network", peptideData, pow=12)
dynamicColors(net1)
}

\author{David L Gibbs}





