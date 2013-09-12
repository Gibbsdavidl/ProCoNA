\name{plotNet}
\alias{plotNet}
\title{plotNet}
\description{Plots the dendrogram and module colors. See ?plotDendroAndColors}
\usage{plotNet(object)}
\arguments{
  \item{object}{The procona network object.}
}

\value{None.}

\examples{
data(ProCoNA_Data)
net1 <- buildProconaNetwork("peptide network", peptideData)
plotNet(net1)
}

\author{David L Gibbs}





