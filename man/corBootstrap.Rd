\name{corBootstrap}
\alias{corBootstrap}
\title{corBootstrap}
\description{This function returns a peptide co-expression network object.}
\usage{corBootstrap(pepdat, networkType = "signed", threshold = 1e-04, 
    tmpSaveFile = T)}
\arguments{
  \item{pepdat}{This variable is the data set with rows as samples and cols as peptides}
  \item{networkType}{Whether the sign is considered in constructing adjacency and TOM}
  \item{threshold}{When to stop running}
  \item{tmpSaveFile}{Should temporary saves be done?}
}

\value{Returns a list of the mean matrix, sd matrix, and the number of resamplings done.}

\author{David L Gibbs}





