\name{adj}
\alias{adj}
\title{adj}
\description{Returns the adjacency matrix.}
\usage{adj(p)}
\arguments{
  \item{p}{The procona network object}
}
\details{
The procona network object contains a number of slots 
which store information relevant to the construction of 
the network. Accessor functions provide direct access to 
the slots. See getSlots("proconaNet") for a complete list.
}
\value{Returns the adjacency matrix.}
\examples{
data(procona)
a <- adj(net)

}
\author{David L Gibbs}





