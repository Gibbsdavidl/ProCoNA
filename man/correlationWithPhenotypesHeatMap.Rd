\name{correlationWithPhenotypesHeatMap}
\alias{correlationWithPhenotypesHeatMap}
\title{correlationWithPhenotypesHeatMap}

\description{Plots a heatmap showing the Pearson correlation of modules with phenotypes.}

\usage{correlationWithPhenotypesHeatMap(net, phenotypes, modules, 
    plotName, title, textSize)}

\arguments{
  \item{net}{The ProCoNA network object.}
  \item{phenotypes}{Matrix of phenotypic traits, can include character strings (converted to factors).}
  \item{modules}{Vector of modules to plot. Default is all modules.}
  \item{plotName}{Name of the saved plot, NULL to show on screen.}
  \item{title}{Plot title.}
  \item{textSize}{The font size of the correlations shown in each module-phenotype pair.}
}

\value{the module eigenvector correlations}

\examples{
data(ProCoNA_Data)
net <- buildProconaNetwork("pepnet", peptideData, pow=12)
n <- length(samples(net))
phenotypes <- matrix(rnorm(10*n), nrow=60)
moduleCors <- correlationWithPhenotypesHeatMap(net, phenotypes, modules = 1:7, 
    plotName = "Phenotype Associations", title = "Module-trait relationships", textSize = 0.5)
}
\author{David L Gibbs}





