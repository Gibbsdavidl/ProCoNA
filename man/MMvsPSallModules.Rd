\name{MMvsPSallModules}
\alias{MMvsPSallModules}
\title{MMvsPSallModules}
\description{Call MMvsPS, producing plots for all modules.}
\usage{MMvsPSallModules(net, peptable, phenoVec, prefixName = "mm_vs_ps_")}
\arguments{
  \item{net}{The procona network object}
  \item{peptable}{The peptide data}
  \item{phenoVec}{The phenotypic trait, as a numeric vector}
  \item{prefixName}{The plot files prefix name.  Writes pdfs.}
}

\value{nothing returned}


\examples{
\dontrun{
# This function outputs a set of pdfs.
data(ProCoNA_Data)
net1 <- buildProconaNetwork("peptide network", peptideData, pow=13)
MMvsPSallModules(net1, peptideData, phenotypes[,5], 1)
}
}

\author{David L Gibbs}





