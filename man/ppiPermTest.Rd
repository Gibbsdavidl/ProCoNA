\name{ppiPermTest}
\alias{ppiPermTest}
\title{ppiPermTest}
\description{Performs a permutation test for enrichment of PPI edges given a database.}
\usage{
ppiPermTest(threshold = 0.33,
iterations = 1000,
pi_edges,
pepinfo,
pepnet,
pepdat,
netRefSeqs)
}
\arguments{
  \item{threshold}{Centrality threshold}
  \item{iterations}{Number of repititions}
  \item{pi_edges}{Must be two columns A-B ... sort out in vivo or in vitro in advance}
  \item{pepinfo}{Maps peptides to proteins ... same format as in ppiTable}
  \item{pepnet}{procona network object}
  \item{pepdat}{the logged data matrix with peptides as columns.}
  \item{netRefSeqs}{the peptides mapped to refseq IDs. named vector.}
}

\value{returns list of test results.}

\author{David L Gibbs}





