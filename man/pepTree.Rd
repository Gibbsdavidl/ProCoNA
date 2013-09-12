\name{pepTree}
\alias{pepTree}
\title{pepTree}
\description{The peptide dendrogram.}
\usage{pepTree(p)}
\arguments{
  \item{p}{the procona obj}
}

\value{Returns the peptide dendrogram, as an hclust object.}

\examples{
\dontrun{
data(ProCoNA_Data)
net1 <- buildProconaNetwork("peptide network", peptideData)
pepTree(net1)
}
}

\author{David L Gibbs}





