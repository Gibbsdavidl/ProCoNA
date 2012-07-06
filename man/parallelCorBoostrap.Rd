\name{parallelCorBoostrap}
\alias{parallelCorBoostrap}
\title{parallelCorBoostrap}
\description{This function returns a peptide co-expression network object.}
\usage{parallelCorBoostrap(pepdat, networkType = "signed", threads = 1, 
    threshold = 1e-04, tmpSaveFiles = T)}
\arguments{
  \item{pepdat}{This variable is the data set with rows as samples and cols as peptides}
  \item{networkType}{Whether the sign is considered in constructing adjacency and TOM}
  \item{threads}{NUMBER OF Threads}
  \item{threshold}{when to stop each thread.}
  \item{tmpSaveFiles}{
}
}

\value{returns the bootstrapped adjacency matrix and sd matrix}

\author{David L Gibbs}





