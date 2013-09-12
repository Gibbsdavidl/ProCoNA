\name{modulePhenotypeCorrelations}
\alias{modulePhenotypeCorrelations}
\title{modulePhenotypeCorrelations}

\description{Computes the relation between the modules and
the phenotypes.}

\usage{modulePhenotypeCorrelations(pnet, phenotypes)}

\arguments{
  \item{pnet}{The peptide net object}
  \item{phenotypes}{The matrix of traits}
}

\value{returns a matrix of correlations between modules and phenotypes.}

\examples{
data(ProCoNA_Data)
net1 <- buildProconaNetwork("peptide network", peptideData, pow=13)
n <- length(samples(net1))
phenotypes <- matrix(rnorm(10*n), nrow=60)
m <- modulePhenotypeCorrelations(net1, phenotypes)

# To plot the heatmap: 
# moduleCors <- correlationWithPhenotypesHeatMap(net1, phenotypes, modules = 1:5, 
#    plot = NULL, title = "Module-trait relationships", textSize = 0.5)

}
\author{David L Gibbs}
